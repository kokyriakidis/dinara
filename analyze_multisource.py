#!/usr/bin/env python3
"""
Analyze multi-source overlap at each locus.

Parses PRE-DEDUP and DELETION CALL lines from dinara logs to answer:
1. How often do multiple sources call the same locus?
2. When they do, how often do they disagree on size?
3. Could a higher-precision source correct a lower-precision source's sizing?

Usage: python3 analyze_multisource.py <results_dir> <truvari_tp_file> <truvari_fn_file>

The truvari files are the tp-call.vcf and fn.vcf from a truvari bench run.
"""

import sys
import os
import re
import glob
from collections import defaultdict

# Source precision from the essentiality analysis (V36i).
SOURCE_PRECISION = {
    "merged-clusters": 0.902,
    "diagonal":        0.807,
    "SA-tag":          0.564,
    "early-CIGAR":     0.486,
    "split-read":      0.472,
    "cluster":         0.325,
    "flank-gap":       0.320,
    "path-based":      0.299,
    "kmer-journey":    0.281,
    "per-read-DEL":    0.281,
    "INV-cluster":     0.194,
    "multi-k":         0.138,
}

def parse_logs(results_dir):
    """Parse all .log files, extract PRE-DEDUP and DELETION CALL lines per locus."""
    loci = {}  # locus_name -> {"pre_dedup": [...], "emitted": [...]}

    for logfile in sorted(glob.glob(os.path.join(results_dir, "*.log"))):
        locus = os.path.basename(logfile).replace(".log", "")
        pre_dedup = []
        emitted = []

        with open(logfile) as f:
            for line in f:
                line = line.strip()
                m = re.match(
                    r'PRE-DEDUP: source=(\S+), size=(\d+), bp=(\d+), reads=(\d+)',
                    line
                )
                if m:
                    pre_dedup.append({
                        "source": m.group(1),
                        "size": int(m.group(2)),
                        "bp": int(m.group(3)),
                        "reads": int(m.group(4)),
                    })
                    continue

                m = re.match(
                    r'>>> DELETION CALL \((\S+)\): size=(\d+)bp, breakpoint=(\d+), reads=(\d+)',
                    line
                )
                if m:
                    emitted.append({
                        "source": m.group(1),
                        "size": int(m.group(2)),
                        "bp": int(m.group(3)),
                        "reads": int(m.group(4)),
                    })

        if pre_dedup or emitted:
            loci[locus] = {"pre_dedup": pre_dedup, "emitted": emitted}

    return loci


def parse_vcf_ids(vcf_file):
    """Extract locus IDs from a VCF file (ID column)."""
    ids = set()
    if not os.path.exists(vcf_file):
        return ids
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                ids.add(fields[2])
    return ids


def size_ratio(a, b):
    return min(a, b) / max(a, b) if max(a, b) > 0 else 0


def analyze(loci, tp_ids, fn_ids):
    # Stats
    total_loci = len(loci)
    multi_source_loci = 0
    single_source_loci = 0
    size_agree = 0
    size_disagree = 0

    # Cases where a high-precision source has a different size than the winner
    consensus_could_help = []
    # Cases where dedup picked a low-precision source over a high-precision one
    wrong_winner = []

    for locus, data in sorted(loci.items()):
        calls = data["pre_dedup"]
        if not calls:
            continue

        sources = set(c["source"] for c in calls)
        if len(sources) <= 1:
            single_source_loci += 1
            continue

        multi_source_loci += 1

        # Group calls by source, pick best (most reads) per source
        best_per_source = {}
        for c in calls:
            src = c["source"]
            if src not in best_per_source or c["reads"] > best_per_source[src]["reads"]:
                best_per_source[src] = c

        # Check size agreement: do all sources agree within ratio 0.9?
        sizes = [c["size"] for c in best_per_source.values()]
        all_agree = True
        for i in range(len(sizes)):
            for j in range(i + 1, len(sizes)):
                if size_ratio(sizes[i], sizes[j]) < 0.9:
                    all_agree = False
                    break
            if not all_agree:
                break

        if all_agree:
            size_agree += 1
        else:
            size_disagree += 1

        # The emitted (post-dedup) winner
        emitted = data["emitted"]
        if not emitted:
            continue

        winner = emitted[0]  # highest read count after dedup
        winner_prec = SOURCE_PRECISION.get(winner["source"], 0)

        # Find highest-precision source that disagrees with winner
        for src, call in sorted(best_per_source.items(),
                                key=lambda x: SOURCE_PRECISION.get(x[0], 0),
                                reverse=True):
            src_prec = SOURCE_PRECISION.get(src, 0)
            if src == winner["source"]:
                continue
            if size_ratio(call["size"], winner["size"]) < 0.9:
                # This high-precision source disagrees on size
                # Determine truth status
                truth_status = "?"
                if locus in tp_ids:
                    truth_status = "TP"
                elif locus in fn_ids:
                    truth_status = "FN"

                consensus_could_help.append({
                    "locus": locus,
                    "winner_source": winner["source"],
                    "winner_size": winner["size"],
                    "winner_reads": winner["reads"],
                    "winner_prec": winner_prec,
                    "alt_source": src,
                    "alt_size": call["size"],
                    "alt_reads": call["reads"],
                    "alt_prec": src_prec,
                    "truth": truth_status,
                })
                break  # only report the highest-precision alternative

    # Print results
    print(f"=== Multi-Source Overlap Analysis ===")
    print(f"Total loci with calls: {total_loci}")
    print(f"Single-source loci:    {single_source_loci}")
    print(f"Multi-source loci:     {multi_source_loci}")
    print(f"  Size agreement:      {size_agree}")
    print(f"  Size disagreement:   {size_disagree}")
    print()

    if consensus_could_help:
        print(f"=== Cases Where Higher-Precision Source Disagrees ({len(consensus_could_help)}) ===")
        print(f"{'Locus':<35} {'Truth':>5} {'Winner':>16} {'WSize':>6} {'WReads':>6} {'WPrec':>5} "
              f"{'AltSrc':>16} {'ASize':>6} {'AReads':>6} {'APrec':>5}")
        print("-" * 140)

        # Count by truth status
        tp_cases = [c for c in consensus_could_help if c["truth"] == "TP"]
        fn_cases = [c for c in consensus_could_help if c["truth"] == "FN"]
        unk_cases = [c for c in consensus_could_help if c["truth"] == "?"]

        for c in sorted(consensus_could_help, key=lambda x: x["truth"]):
            print(f"{c['locus']:<35} {c['truth']:>5} {c['winner_source']:>16} "
                  f"{c['winner_size']:>6} {c['winner_reads']:>6} {c['winner_prec']:>5.2f} "
                  f"{c['alt_source']:>16} {c['alt_size']:>6} {c['alt_reads']:>6} {c['alt_prec']:>5.2f}")

        print()
        print(f"Summary: {len(tp_cases)} TP, {len(fn_cases)} FN, {len(unk_cases)} unknown")
        print(f"  FN cases = wrong-size winner where a better source had a different size")
        print(f"  These are the cases where consensus could potentially fix the call")

    # Source co-occurrence matrix
    print()
    print("=== Source Co-occurrence (multi-source loci) ===")
    cooccur = defaultdict(int)
    for locus, data in loci.items():
        sources = set(c["source"] for c in data["pre_dedup"])
        if len(sources) <= 1:
            continue
        for s1 in sources:
            for s2 in sources:
                if s1 <= s2:
                    cooccur[(s1, s2)] += 1

    all_sources = sorted(set(s for pair in cooccur for s in pair))
    print(f"{'':>16}", end="")
    for s in all_sources:
        print(f" {s[:8]:>8}", end="")
    print()
    for s1 in all_sources:
        print(f"{s1:>16}", end="")
        for s2 in all_sources:
            key = (min(s1, s2), max(s1, s2))
            print(f" {cooccur.get(key, 0):>8}", end="")
        print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results_dir> [tp-call.vcf] [fn.vcf]")
        sys.exit(1)

    results_dir = sys.argv[1]
    tp_file = sys.argv[2] if len(sys.argv) > 2 else ""
    fn_file = sys.argv[3] if len(sys.argv) > 3 else ""

    tp_ids = parse_vcf_ids(tp_file) if tp_file else set()
    fn_ids = parse_vcf_ids(fn_file) if fn_file else set()

    loci = parse_logs(results_dir)
    analyze(loci, tp_ids, fn_ids)
