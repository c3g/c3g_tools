#!/usr/bin/env python3
"""
s3_md5_metadata.py
==================
Compute MD5 of S3 objects via streaming (no full download to disk/RAM)
and update each object's metadata in-place using a server-side copy.

Architecture
------------
1. List all objects upfront (fast: ~1s for thousands of objects).
   This gives full workload visibility before any processing starts.

2. Split into two lists:
     Fast lane  -- small objects (< --large-object-mb)
                   multiprocessing.Pool, many workers, small chunks.
                   Optimised for high object/s throughput.
     Slow lane  -- large objects (>= --large-object-mb)
                   ThreadPoolExecutor, fewer workers, larger chunks.
                   Sorted largest-first so the longest tasks start
                   immediately, minimising straggler idle time at the end.

3. Both lanes run concurrently in separate threads against their
   pre-built task lists. No queues, no events, no sharding complexity.

Efficiency notes
----------------
- RAM   : full object list in memory (O(n) key strings, typically small).
          Chunk buffers: workers x chunk-mb per lane.
- CPU   : hashlib uses OpenSSL's C MD5 implementation.
          Fast lane uses multiprocessing to bypass GIL for hashing.
          Slow lane uses threading — workers spend >95% of time in network
          I/O where the GIL is released, so threads are as fast as processes.
- I/O   : no disk writes.
- BW    : one GET per object for hashing; metadata update is a server-side
          CopyObject (no content re-upload). Single HEAD per object reused
          for skip check and copy step.

TCP connection reuse
--------------------
All clients use tcp_keepalive=True and pool sizes matched to actual
concurrency to prevent urllib3 "pool is full" connection discard.

Credentials
-----------
boto3 default chain: env vars -> ~/.aws/credentials -> IAM role -> SSO.

Usage
-----
# Basic
python s3_md5_metadata.py --bucket my-bucket --endpoint-url https://...

# Tuned for Ceph RGW on Juno
python s3_md5_metadata.py --bucket my-bucket \\
    --endpoint-url https://objets.juno.calculquebec.ca \\
    --fast-workers 6  --fast-chunk-mb 16 \\
    --slow-workers 8  --slow-chunk-mb 64 \\
    --large-object-mb 500

# Force recompute, verbose
python s3_md5_metadata.py --bucket my-bucket --force --progress-every 50
"""

import argparse
import hashlib
import logging
import multiprocessing
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
from botocore.config import Config

log = logging.getLogger(__name__)

METADATA_KEY        = "md5chksum"
COPY_OBJECT_LIMIT   = 5 * 1024**3    # 5 GiB — single CopyObject limit
MULTIPART_PART_SIZE = 64 * 1024**2   # 64 MiB per multipart part


# ---------------------------------------------------------------------------
# boto3 client factory
# ---------------------------------------------------------------------------

def _make_client(endpoint_url: str | None, pool_connections: int = 4):
    """
    pool_connections must match the actual concurrency of the caller to
    avoid urllib3 discarding connections ("pool is full" warning).
    tcp_keepalive keeps connections alive between tasks.
    """
    cfg = Config(
        max_pool_connections=pool_connections,
        tcp_keepalive=True,
        retries={"max_attempts": 5, "mode": "adaptive"},
    )
    return boto3.client("s3", endpoint_url=endpoint_url, config=cfg)


# ---------------------------------------------------------------------------
# Listing — returns two sorted lists
# ---------------------------------------------------------------------------

def list_objects(
    endpoint_url: str | None,
    bucket: str,
    prefix: str,
    threshold_bytes: int,
    skip_existing: bool,
) -> tuple[list, list]:
    """
    Page through the bucket and return:
      fast_tasks -- [(key, size), ...] for objects below threshold, unsorted
      slow_tasks -- [(key, size), ...] for objects at/above threshold,
                    sorted largest-first so straggler tasks start immediately.

    When skip_existing=True, objects that already carry the MD5 metadata key
    are excluded here — saving a HEAD call later for each skipped object.
    ListObjectsV2 does not return user metadata, so we cannot skip at list
    time; the HEAD is still needed per object. We keep this simple: include
    all objects and let the worker do the skip check.
    """
    client    = _make_client(endpoint_url, pool_connections=2)
    paginator = client.get_paginator("list_objects_v2")
    fast, slow = [], []

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            item = (obj["Key"], obj["Size"])
            if obj["Size"] >= threshold_bytes:
                slow.append(item)
            else:
                fast.append(item)

    # Largest objects first — workers pick up the heaviest tasks immediately,
    # minimising idle time at the tail of the slow lane.
    slow.sort(key=lambda x: x[1], reverse=True)

    return fast, slow


# ===========================================================================
# FAST LANE — multiprocessing
# ===========================================================================

_fast_s3     = None
_fast_chunk  = None
_fast_skip   = None
_fast_bucket = None


def _fast_worker_init(bucket: str, chunk_size: int, skip_existing: bool,
                      endpoint_url: str | None):
    global _fast_s3, _fast_chunk, _fast_skip, _fast_bucket
    _fast_s3     = _make_client(endpoint_url, pool_connections=2)
    _fast_chunk  = chunk_size
    _fast_skip   = skip_existing
    _fast_bucket = bucket


def _fast_process(task: tuple) -> dict:
    key, size = task
    return _do_object(key, size, _fast_s3, _fast_chunk, _fast_skip, _fast_bucket)


# ===========================================================================
# SLOW LANE — threading
# ===========================================================================

_slow_local = threading.local()


def _get_slow_client(endpoint_url: str | None):
    if not hasattr(_slow_local, "s3"):
        _slow_local.s3 = _make_client(endpoint_url, pool_connections=2)
    return _slow_local.s3


def _make_slow_worker(bucket: str, chunk_size: int, skip_existing: bool,
                      endpoint_url: str | None):
    def _worker(task: tuple) -> dict:
        key, size = task
        return _do_object(key, size,
                          _get_slow_client(endpoint_url),
                          chunk_size, skip_existing, bucket)
    return _worker


# ===========================================================================
# Shared object processing
# ===========================================================================

def _do_object(key: str, size: int, s3, chunk_size: int,
               skip_existing: bool, bucket: str) -> dict:
    result = {"key": key, "status": None, "md5": None, "size": size}
    try:
        head = s3.head_object(Bucket=bucket, Key=key)

        if skip_existing and METADATA_KEY in head.get("Metadata", {}):
            result["status"] = "skipped"
            result["md5"]    = head["Metadata"][METADATA_KEY]
            log.debug("SKIP  %s", key)
            return result

        md5_hex = _stream_md5(key, s3, chunk_size, bucket)
        _update_metadata(key, md5_hex, head, s3, bucket)

        result["status"] = "updated"
        result["md5"]    = md5_hex
        log.debug("OK    %s  md5=%s", key, md5_hex)

    except Exception as exc:  # noqa: BLE001
        result["status"] = "error"
        result["error"]  = str(exc)
        log.error("FAIL  %s  --  %s", key, exc)

    return result


def _stream_md5(key: str, s3, chunk_size: int, bucket: str) -> str:
    hasher   = hashlib.md5()
    response = s3.get_object(Bucket=bucket, Key=key)
    for chunk in response["Body"].iter_chunks(chunk_size=chunk_size):
        hasher.update(chunk)
    return hasher.hexdigest()


def _update_metadata(key: str, md5_hex: str, head: dict, s3, bucket: str) -> None:
    meta = dict(head.get("Metadata", {}))
    meta[METADATA_KEY] = md5_hex
    extra = {
        "ContentType": head.get("ContentType", "binary/octet-stream"),
        "StorageClass": head.get("StorageClass", "STANDARD"),
    }
    if head.get("ContentLength", 0) > COPY_OBJECT_LIMIT:
        _multipart_copy(key, meta, extra, head["ContentLength"], s3, bucket)
    else:
        s3.copy_object(
            Bucket=bucket, Key=key,
            CopySource={"Bucket": bucket, "Key": key},
            Metadata=meta, MetadataDirective="REPLACE",
            **extra,
        )


def _multipart_copy(key: str, meta: dict, extra: dict,
                    content_length: int, s3, bucket: str) -> None:
    upload_id = None
    try:
        resp      = s3.create_multipart_upload(
            Bucket=bucket, Key=key, Metadata=meta, **extra)
        upload_id = resp["UploadId"]
        parts     = []
        for part_num, start in enumerate(
                range(0, content_length, MULTIPART_PART_SIZE), 1):
            end  = min(content_length, start + MULTIPART_PART_SIZE) - 1
            part = s3.upload_part_copy(
                Bucket=bucket, Key=key,
                PartNumber=part_num, UploadId=upload_id,
                CopySource={"Bucket": bucket, "Key": key},
                CopySourceRange=f"bytes={start}-{end}",
            )
            parts.append({"PartNumber": part_num,
                          "ETag": part["CopyPartResult"]["ETag"]})
        s3.complete_multipart_upload(
            Bucket=bucket, Key=key,
            UploadId=upload_id,
            MultipartUpload={"Parts": parts},
        )
    except Exception:
        if upload_id is not None:
            s3.abort_multipart_upload(Bucket=bucket, Key=key, UploadId=upload_id)
        raise


# ===========================================================================
# Lane stats (thread-safe)
# ===========================================================================

class LaneStats:
    def __init__(self, name: str, total_objects: int,
                 progress_every: int, t_start: float):
        self.name           = name
        self.total_objects  = total_objects   # known upfront from listing
        self.progress_every = progress_every
        self.t_start        = t_start
        self._lock          = threading.Lock()
        self.done           = 0
        self.bytes_done     = 0
        self.updated        = 0
        self.skipped        = 0
        self.error          = 0
        self.t_first        = None
        self.t_end          = None

    def record(self, result: dict) -> None:
        with self._lock:
            if self.t_first is None:
                self.t_first = time.monotonic()
            self.done       += 1
            self.bytes_done += result.get("size", 0)
            self.updated    += result["status"] == "updated"
            self.skipped    += result["status"] == "skipped"
            self.error      += result["status"] == "error"
            should_log = (self.done % self.progress_every == 0)
        if should_log:
            self._log_progress()

    def finish(self) -> None:
        with self._lock:
            if self.t_end is None:
                self.t_end = time.monotonic()
        self._log_summary()

    def _elapsed(self) -> float:
        end = self.t_end or time.monotonic()
        return max(0.001, end - (self.t_first or self.t_start))

    def _log_progress(self) -> None:
        e    = self._elapsed()
        mibs = self.bytes_done / e / 1024 / 1024
        pct  = 100 * self.done / self.total_objects if self.total_objects else 0
        log.info("[%s] %d/%d (%.0f%%)  |  %.1f obj/s  |  %.1f MiB/s  |  "
                 "updated=%d  skipped=%d  errors=%d",
                 self.name, self.done, self.total_objects, pct,
                 self.done / e, mibs,
                 self.updated, self.skipped, self.error)

    def _log_summary(self) -> None:
        e    = self._elapsed()
        mibs = self.bytes_done / e / 1024 / 1024
        log.info("[%s] Done -- %d objects in %.1fs  |  %.1f obj/s  |  "
                 "%.1f MiB/s (%.2f GiB/s)  |  "
                 "updated=%d  skipped=%d  errors=%d",
                 self.name, self.done, e, self.done / e,
                 mibs, mibs / 1024,
                 self.updated, self.skipped, self.error)


# ===========================================================================
# CLI
# ===========================================================================

def parse_args():
    cpu_count    = multiprocessing.cpu_count()
    fast_default = max(1, cpu_count // 2)
    slow_default = max(1, cpu_count // 8)

    parser = argparse.ArgumentParser(
        description="Stream-hash S3 objects and store MD5 in object metadata.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bucket",       required=True)
    parser.add_argument("--prefix",       default="",
                        help="Key prefix filter ('' = whole bucket).")
    parser.add_argument("--endpoint-url", default=None,
                        help="Custom S3 endpoint URL (e.g. Ceph RGW).")
    parser.add_argument("--force",        action="store_true",
                        help="Recompute MD5 even if metadata key already exists.")
    parser.add_argument("--progress-every", type=int, default=100,
                        help="Log progress every N objects per lane.")
    parser.add_argument("--log-level",    default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    lane = parser.add_argument_group("lane tuning")
    lane.add_argument("--large-object-mb", type=float, default=500.0,
                      help="Size threshold MiB. Objects >= this go to slow lane.")
    lane.add_argument("--fast-workers",    type=int,   default=fast_default,
                      help="Worker processes for fast lane.")
    lane.add_argument("--fast-chunk-mb",   type=float, default=16.0,
                      help="Streaming chunk size MiB for fast lane.")
    lane.add_argument("--slow-workers",    type=int,   default=slow_default,
                      help="Worker threads for slow lane.")
    lane.add_argument("--slow-chunk-mb",   type=float, default=64.0,
                      help="Streaming chunk size MiB for slow lane.")
    return parser.parse_args()


# ===========================================================================
# Main
# ===========================================================================

def main():
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    threshold_bytes = int(args.large_object_mb * 1024 * 1024)
    skip_existing   = not args.force
    fast_chunk      = int(args.fast_chunk_mb * 1024 * 1024)
    slow_chunk      = int(args.slow_chunk_mb * 1024 * 1024)
    peak_ram        = (args.fast_workers * args.fast_chunk_mb
                       + args.slow_workers * args.slow_chunk_mb)

    log.info(
        "Fast lane: %d workers x %.0f MiB chunks  |  "
        "Slow lane: %d threads x %.0f MiB chunks  |  "
        "Threshold: %.0f MiB  |  Peak RAM approx %.0f MiB",
        args.fast_workers, args.fast_chunk_mb,
        args.slow_workers, args.slow_chunk_mb,
        args.large_object_mb, peak_ram,
    )
    log.info("Listing s3://%s/%s ...", args.bucket, args.prefix)

    # Connectivity check.
    probe = _make_client(args.endpoint_url)
    try:
        probe.head_bucket(Bucket=args.bucket)
    except Exception:
        log.exception("Failed to connect to bucket '%s'.", args.bucket)
        raise

    t_list_start = time.monotonic()
    fast_tasks, slow_tasks = list_objects(
        args.endpoint_url, args.bucket, args.prefix,
        threshold_bytes, skip_existing,
    )
    t_list = time.monotonic() - t_list_start

    total = len(fast_tasks) + len(slow_tasks)
    log.info(
        "Listed %d objects in %.1fs -- "
        "fast lane: %d  |  slow lane: %d (largest-first)",
        total, t_list, len(fast_tasks), len(slow_tasks),
    )

    if total == 0:
        log.info("Nothing to do.")
        return 0

    t_start       = time.monotonic()
    combined      = {"total": 0, "updated": 0, "skipped": 0, "error": 0}
    combined_lock = threading.Lock()

    slow_progress = max(1, args.progress_every // 10)
    fast_stats    = LaneStats("fast", len(fast_tasks), args.progress_every, t_start)
    slow_stats    = LaneStats("slow", len(slow_tasks), slow_progress,       t_start)

    def run_fast_lane(pool):
        if not fast_tasks:
            fast_stats.finish()
            return
        for result in pool.imap_unordered(
            _fast_process, fast_tasks,
            chunksize=max(1, args.fast_workers * 4),
        ):
            fast_stats.record(result)
            with combined_lock:
                combined[result["status"]] += 1
                combined["total"]          += 1
        fast_stats.finish()

    def run_slow_lane(executor):
        if not slow_tasks:
            slow_stats.finish()
            return
        slow_worker = _make_slow_worker(
            args.bucket, slow_chunk, skip_existing, args.endpoint_url)
        # Submit all tasks — executor enforces max_workers concurrency.
        # Tasks are already sorted largest-first so the pool picks them up
        # in the right order.
        futures = {executor.submit(slow_worker, task): task
                   for task in slow_tasks}
        for fut in as_completed(futures):
            result = fut.result()
            slow_stats.record(result)
            with combined_lock:
                combined[result["status"]] += 1
                combined["total"]          += 1
        slow_stats.finish()

    try:
        with (
            multiprocessing.Pool(
                processes=args.fast_workers,
                initializer=_fast_worker_init,
                initargs=(args.bucket, fast_chunk, skip_existing, args.endpoint_url),
            ) as fast_pool,
            ThreadPoolExecutor(max_workers=args.slow_workers) as slow_executor,
        ):
            fast_thread = threading.Thread(
                target=run_fast_lane, args=(fast_pool,), name="fast-lane")
            slow_thread = threading.Thread(
                target=run_slow_lane, args=(slow_executor,), name="slow-lane")

            fast_thread.start()
            slow_thread.start()
            fast_thread.join()
            slow_thread.join()

    except KeyboardInterrupt:
        elapsed = time.monotonic() - t_start
        log.warning(
            "Interrupted after %d objects in %.1fs  |  "
            "updated=%d  skipped=%d  errors=%d",
            combined["total"], elapsed,
            combined["updated"], combined["skipped"], combined["error"],
        )
        return 130

    elapsed = time.monotonic() - t_start
    rate    = combined["total"] / elapsed if elapsed > 0 else 0
    log.info(
        "Total -- %d objects in %.1fs (%.1f obj/s)  |  "
        "updated=%d  skipped=%d  errors=%d",
        combined["total"], elapsed, rate,
        combined["updated"], combined["skipped"], combined["error"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
