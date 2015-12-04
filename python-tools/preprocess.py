"""
This will standardize input file genotype fields to have AD=REF,ALT1,ALT2 ... (a
la GTK).

For platypus:  --depth NR --alt-depth NV
    varscan:   --ref-depth RD --alt-depth AD
    samtools:  --ref-depth DP4 (untested)

"""
import argparse
import sys
import gzip
import re

def xopen(f):
    if f.endswith(".gz"): return gzip.open(f)
    else: return open(f)

def dp4(sample, dp4_i):

    part = sample[dp4_i].split(",")

    depths = [int(part[0]) + int(part[1])]
    depths.append(int(part[2]) + int(part[3]))
    sample.append(",".join(str(d) for d in depths))


def set_depths(sample, ref_i, alt_i, dp_i, add_dp, fmt=None):

    if sample[0] in ("./.", "."):
        # if ungenotyped, just truncate the sample filed to be "./."
        sample[:] = sample[:1]
        return

    depths = []
    if ref_i:
        depths.append(sample[ref_i])
    else:
        depths.append(None)
    if alt_i:
        depths.append(sample[alt_i])

    dp = sample[dp_i]
    if "," in dp:
        #assert len(set(dp.split(","))) == 1, (dp, sample)
        dp = map(int, dp.split(",")) 

        sample[dp_i] = dp

    """
    GT:GL:GOF:GQ:NR:NV      1/0:-1,-1,-1:13:3:80,76:37,1
    $ ##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
    $ ##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
    """
    # only had total depth and alt_depth. infer ref depth, a la platypus
    if not ref_i:
        ref_depth = int(dp) - max(int(x) for x in depths[1].split(","))
        assert ref_depth >= 0, (depths, dp)
        depths[0] = str(ref_depth)

    assert not "," in dp

    assert len(depths) > 1, (depths, sample)
    sample.append(",".join(depths))

    # if the caller doesn't use "DP", we add it here.
    if add_dp:
        sample.append(sample[dp_i])


def preprocess(vcf, ref_depth, alt_depth, depth="DP"):

    n_samples = 0
    for line in xopen(vcf):

        if line.startswith("#"):
            if line.startswith("##FORMAT"):
                if "ID=%s," % ref_depth in line or ref_depth is None \
                        and "ID=%s," % alt_depth in line:
                    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">"
                
                if depth != "DP" and "ID=%s" % depth in line:
                    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"read depths\">"
            if line.startswith("##FORMAT=<ID=AD,"):
                line = line.replace("##FORMAT=<ID=AD,", "##FORMAT=<ID=AD_O,")
            print line,

            if line.startswith("#CHROM"):
                n_samples = len(line.split("\t")) - 9

        else: # past VCF header.
            toks = line.rstrip().rsplit("\t", n_samples + 1)
            toks, fmt, samples = toks[0], toks[1].split(":"), [x.split(":") for x in toks[2:]]

            ref_i = fmt.index(ref_depth) if ref_depth is not None else None
            alt_i = fmt.index(alt_depth) if alt_depth is not None else None
            dp_i = fmt.index(depth)
            
            if "AD" in fmt:
                fmt = [x if x != "AD" else "AD_O" for x in fmt]

            if "DP4" == ref_depth:
                [dp4(s, ref_i) for s in samples]
            else:
                [set_depths(s, ref_i, alt_i, dp_i, depth != "DP") for s in samples]
            fmt.append("AD")

            if depth != "DP": fmt.append("DP")
    
            #print samples[0]
            #print len(samples[0])
            #print len(fmt)
            #print samples[0][0]
            #assert len(samples[0]) == len(fmt) or samples[0][0] == "."

            print "\t".join([toks, ":".join(fmt), "\t".join(":".join(s) for s in samples)])

def main(args):

    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--ref-depth", help="field containing ref-depths")
    p.add_argument("--alt-depth", help="field containing alt-depths (if this is "
            " empty, ref-depth is assumed to contain all depths with ref first)")
    p.add_argument("--depth", help="field containing total depths", default="DP")
    p.add_argument("vcf", help="vcf to normalize")

    a = p.parse_args(args)

    preprocess(a.vcf, a.ref_depth, a.alt_depth, a.depth)


if __name__ == "__main__":

    main(sys.argv[1:])
