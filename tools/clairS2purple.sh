#!/bin/env bash

# Convert ClairS somatic VCF to PURPLE-compatible format
# Re-implemented as a C3G_TOOLS script

set -e
set -o pipefail

usage() {
    echo "Usage: $(basename $0) -i <input.vcf.gz> -t <TUMOR> -n <NORMAL>" >&2
    echo "  -i  Input ClairS somatic VCF (compressed)" >&2
    echo "  -t  Tumor sample name" >&2
    echo "  -n  Normal sample name" >&2
    exit 1
}

INPUT_VCF=""
OUTPUT_VCF=""
TUMOR_SAMPLE=""
NORMAL_SAMPLE=""

while getopts ":i:t:n:" opt; do
    case $opt in
        i) INPUT_VCF=$OPTARG ;;
        t) TUMOR_SAMPLE=$OPTARG ;;
        n) NORMAL_SAMPLE=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_VCF" ]] || [[ -z "$TUMOR_SAMPLE" ]] || [[ -z "$NORMAL_SAMPLE" ]]; then
    echo "ERROR: Missing required arguments" >&2
    usage
fi

if [[ ! -s "$INPUT_VCF" ]]; then
    echo "ERROR: Input VCF file does not exist or is empty: $INPUT_VCF" >&2
    exit 1
fi

echo "Converting ClairS VCF to PURPLE format..." >&2

# Compress and convert VCF: reconstruct allele depths for SNVs and indels
zcat "$INPUT_VCF" | awk -v tumor="$TUMOR_SAMPLE" -v normal="$NORMAL_SAMPLE" '
BEGIN { FS = OFS = "\t" }

# Preserve VCF metadata
/^##/ { print; next }

# Update header with new FORMAT fields
/^#CHROM/ {
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">"
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">"
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, tumor, normal
    next
}

# Process variants: reconstruct allele depths for SNVs and indels
{
    REF = $4; ALT = $5
    TR = TA = NR = NA = 0

    split($9, FMT, ":")
    split($10, VAL, ":")
   
   # clears the lookup table for each new variant line, so values from the previous VCF row do not leak into the next one.
    delete DATA
    for (i = 1; i <= length(FMT); i++) DATA[FMT[i]] = VAL[i]

    if (length(REF) == 1 && length(ALT) == 1) {
        # SNV: reconstruct from nucleotide counts (ClairS AD field is incorrect)
        AU = DATA["AU"] + 0; CU = DATA["CU"] + 0; GU = DATA["GU"] + 0; TU = DATA["TU"] + 0
        NAU = DATA["NAU"] + 0; NCU = DATA["NCU"] + 0; NGU = DATA["NGU"] + 0; NTU = DATA["NTU"] + 0

        if (REF == "A") { TR = AU; NR = NAU } else if (REF == "C") { TR = CU; NR = NCU } else if (REF == "G") { TR = GU; NR = NGU } else if (REF == "T") { TR = TU; NR = NTU }
        if (ALT == "A") { TA = AU; NA = NAU } else if (ALT == "C") { TA = CU; NA = NCU } else if (ALT == "G") { TA = GU; NA = NGU } else if (ALT == "T") { TA = TU; NA = NTU }

        TDP = TR + TA; NDP = NR + NA
    } else {
        # Indel: use AD field (more reliable for indels)
        TDP = DATA["DP"] + 0; NDP = DATA["NDP"] + 0

        split(DATA["AD"], TAD, ","); TA = TAD[2] + 0; TR = TDP - TA
        split(DATA["NAD"], NAD, ","); NA = NAD[2] + 0; NR = NDP - NA
    }

    # Determine normal genotype based on ALT allele fraction (threshold 0.2)
    GT = DATA["GT"]
    if (NDP > 0) AF = NA / NDP; else AF = 0
    NORM_GT = AF > 0.2 ? "0/1" : "0/0"

    TUMOR_OUT = GT ":" TR "," TA ":" TDP
    NORMAL_OUT = NORM_GT ":" NR "," NA ":" NDP
    print $1, $2, $3, $4, $5, $6, $7, $8, "GT:AD:DP", TUMOR_OUT, NORMAL_OUT
}'

echo "Done. Output: $OUTPUT_VCF" >&2
