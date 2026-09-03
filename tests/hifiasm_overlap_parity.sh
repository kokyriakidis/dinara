#!/usr/bin/env bash
# Overlap-generation parity check: does the vendored hifiasmCandidates fork
# (external/hifiasmCandidates, what dinara actually links against) compute the
# same candidate overlaps as real upstream hifiasm for the same reads and
# parameters?
#
# The fork intentionally runs NO error correction (see
# external/hifiasmCandidates/candidates.cpp), so it is not meaningful to diff
# against real hifiasm's normal (corrected-read) *.ovlp.paf output -- reads
# genuinely differ after correction. Real hifiasm's only raw-read,
# single-pass overlap mode is --dbg-ovec, which unconditionally skips the
# k-mer filter table regardless of -f/-F (see ha_ec_dbg() in Assembly.cpp).
# So this script runs two checks instead of one:
#
#   Test A: fork, dinara's real settings (filter ON, default -f/-D/--max-kocc)
#           vs. fork, filter OFF (-F)
#           -- regression guard: catches accidental filter misconfiguration.
#           Expected to match whenever the input has no k-mers extreme enough
#           to trip the high-occurrence filter; if it legitimately diverges
#           for a larger/more repetitive input, that's not itself a fork bug,
#           just means Test B (below) is the one that actually matters there.
#
#   Test B: fork, filter OFF (-F) vs. real hifiasm, --dbg-ovec (always no
#           filter)
#           -- the real parity check: same raw reads, same k/w, no-HPC,
#           no filter, single pass, base-level alignment. Any mismatch here
#           is a genuine divergence between the fork and upstream hifiasm's
#           actual overlap algorithm.
#
# Usage:
#   tests/hifiasm_overlap_parity.sh <reads.fastq> [k] [w] [threads]
#
# Env overrides:
#   FORK_HIFIASM  path to the built external/hifiasmCandidates/hifiasm CLI
#                 (default: <repo>/external/hifiasmCandidates/hifiasm)
#   REAL_HIFIASM  path to a real upstream hifiasm build
#                 (default: /home/kokyriakidis/Downloads/hifiasm/hifiasm)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

READS="${1:?usage: $0 <reads.fastq> [k] [w] [threads]}"
K="${2:-50}"
W="${3:-50}"
THREADS="${4:-4}"

FORK_HIFIASM="${FORK_HIFIASM:-$REPO_ROOT/external/hifiasmCandidates/hifiasm}"
REAL_HIFIASM="${REAL_HIFIASM:-/home/kokyriakidis/Downloads/hifiasm/hifiasm}"

for exe in "$FORK_HIFIASM" "$REAL_HIFIASM"; do
    if [[ ! -x "$exe" ]]; then
        echo "ERROR: not an executable: $exe" >&2
        exit 1
    fi
done

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

echo "[parity] reads=$READS k=$K w=$W threads=$THREADS"
echo "[parity] fork=$FORK_HIFIASM"
echo "[parity] real=$REAL_HIFIASM"

# Keyed as: qname qlen qs qe strand tname tlen ts te <cigar, =/X collapsed to M>
# so fork's extended CIGAR (=/X) and real's basic CIGAR (M) compare equal, and
# the (nmatch, alnlen) columns -- present in the fork's PAF but absent from
# real's debug PAF -- are dropped from the comparison key entirely.
key_fork() {
    awk 'BEGIN{OFS="\t"}{
        cg=$13; sub(/^cg:Z:/,"",cg); gsub(/[=X]/,"M",cg);
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,cg
    }' "$1" | sort
}
key_real() {
    awk 'BEGIN{OFS="\t"}{
        cg=$11; sub(/^cg:Z:/,"",cg);
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,cg
    }' "$1" | sort
}

report() {
    local label="$1" a="$2" b="$3"
    local only_a only_b n_a n_b
    n_a=$(wc -l < "$a"); n_b=$(wc -l < "$b")
    only_a=$(comm -23 "$a" "$b" | wc -l)
    only_b=$(comm -13 "$a" "$b" | wc -l)
    echo "[$label] left=$n_a right=$n_b only_left=$only_a only_right=$only_b"
    if [[ "$only_a" -eq 0 && "$only_b" -eq 0 ]]; then
        echo "[$label] PASS: exact match"
        return 0
    else
        echo "[$label] FAIL: mismatch -- sample only_left:"
        comm -23 "$a" "$b" | head -5
        echo "[$label] FAIL: mismatch -- sample only_right:"
        comm -13 "$a" "$b" | head -5
        return 1
    fi
}

status=0

echo
echo "=== Test A: fork default filter vs fork -F (no filter) ==="
"$FORK_HIFIASM" -k "$K" -w "$W" -t "$THREADS" -o a_filtered "$READS" \
    > a_filtered.log 2>&1
"$FORK_HIFIASM" -k "$K" -w "$W" -F -t "$THREADS" -o a_nofilter "$READS" \
    > a_nofilter.log 2>&1
key_fork a_filtered.ovlp.paf > a_filtered.key
key_fork a_nofilter.ovlp.paf > a_nofilter.key
report "Test A" a_filtered.key a_nofilter.key || status=1

echo
echo "=== Test B: fork -F (no filter) vs real hifiasm --dbg-ovec ==="
"$REAL_HIFIASM" -k "$K" -w "$W" -t "$THREADS" --dbg-ovec -o b_real "$READS" \
    > b_real.log 2>&1
key_fork a_nofilter.ovlp.paf > b_fork.key
key_real b_real.ovlp.paf > b_real.key
report "Test B" b_fork.key b_real.key || status=1

echo
if [[ "$status" -eq 0 ]]; then
    echo "[parity] ALL TESTS PASSED: fork overlaps match real hifiasm exactly."
else
    echo "[parity] FAILED: see above. Logs kept for inspection would normally" \
         "be cleaned up by the trap -- rerun with 'set +e; trap - EXIT' at the" \
         "top if you need to inspect $WORKDIR contents."
fi
exit "$status"
