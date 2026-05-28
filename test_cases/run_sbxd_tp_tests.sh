#!/bin/bash
# Run dinara SV detection on all 80 SBX-D TP test cases.
# Usage: bash test_cases/run_sbxd_tp_tests.sh [dinara_binary]
#
# 80 cases: 10 per size bin × 4 DEL bins × 4 INS bins
# Source: SBX-D.30X.bam (Roche 2x250bp, GRCh38, ~30x coverage)

DINARA=${1:-build/Executable/dinara}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CASES_DIR="$SCRIPT_DIR/sbxd_tp"

if [ ! -x "$DINARA" ]; then
    echo "Error: dinara binary not found at $DINARA"
    echo "Usage: $0 [path/to/dinara]"
    exit 1
fi

echo "=== SBX-D TP TEST CASES (80 cases) ==="
echo "Binary: $DINARA"

for bin in DEL_small_lt100bp DEL_medium_100_500bp DEL_large_500_1000bp DEL_xlarge_gt1000bp \
           INS_small_lt100bp INS_medium_100_500bp INS_large_500_1000bp INS_xlarge_gt1000bp; do
    echo ""
    echo "======== $bin ========"
    for dir in "$CASES_DIR/$bin"/chr*; do
        name=$(basename "$dir")
        true_type=$(grep "SV Type" "$dir/info.txt" | awk '{print $3}')
        true_size=$(grep "SV Length" "$dir/info.txt" | awk '{print $3}')
        echo ""
        echo "--- $name (true: $true_type $true_size) ---"
        outdir=$(mktemp -d)
        "$DINARA" --command svanchors \
            --reference "$dir/reference.fa" \
            --input "$dir/reads.fa" \
            --assemblyDirectory "$outdir" \
            --Reads.minReadLength 50 --Kmers.k 10 --Kmers.minimizerW 6 \
            --bam "$dir/region.bam" 2>&1 \
            | grep -E ">>> |SA-tag refine|Rescued|INS suppressed"
        rm -rf "$outdir"
    done
done
