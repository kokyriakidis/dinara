#!/usr/bin/env python3
# Gate 1 diagnostic: measure collapse/tangle in a Shasta2AnchorGraph GFA.
#
# Reads a GFA written by Shasta2AnchorGraph::writeGfa:
#   S <vid> * LN:i:1 wn:i:<window> ws:Z:<fw|rc>
#   L <src> + <dst> + 0M RC:i:<cov> [sp:i: sn:i: sr:i:]
#
# Reports baseline metrics for the de-novo detangling work:
#   - vertices/edges, distinct normalized windows
#   - vertices-per-window distribution (collapse proxy)
#   - intra- vs inter-window edges
#   - per-window predecessor/successor counts (tangle windows)
#   - strand-strand edges (fw<->rc of same window: hairpin/collapse signature)
#   - hairpin windows (connect to their own normalized window)
#   - span-product (sp*sn) distribution + log2 histogram (drives T threshold)
#   - shared-read (sr) distribution (drives N threshold)
#
# Usage: python3 Gate1AnchorGraphCollapseStats.py Shasta2AnchorGraph.gfa [edges.csv]
#   If edges.csv is given, a per-edge dump of sp,sn,product,sr is written there
#   so T can be re-derived later without re-running.

import re
import sys


def pct(sorted_vals, q):
    n = len(sorted_vals)
    if n == 0:
        return 0
    i = int(q * (n - 1) + 0.5)
    i = max(0, min(i, n - 1))
    return sorted_vals[i]


def log2bucket(v):
    if v <= 0:
        return "0"
    lo = 1
    while lo * 2 <= v:
        lo *= 2
    hi = lo * 2
    return "%d-%d" % (lo, hi - 1)  # e.g. 1-1, 2-3, 4-7, 8-15


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: %s <anchorgraph.gfa> [edges.csv]" % sys.argv[0])
    gfa_path = sys.argv[1]
    csv_path = sys.argv[2] if len(sys.argv) > 2 else None

    v_window = {}   # vid -> normalized window
    v_strand = {}   # vid -> fw|rc
    win_verts = {}  # window -> count of vertices (any strand)
    s_lines = 0
    l_lines = 0

    # Edge bookkeeping at normalized-window granularity.
    succ = {}  # win -> {succWin: count}  (different window only)
    pred = {}  # win -> {predWin: count}
    intra = 0
    inter_same_win_diff_strand = 0
    inter_diff_win = 0

    edges = []          # (srcWin, srcStrand, dstWin, dstStrand)
    span_products = []   # sp*sn per inter-window edge that carries sp/sn
    shared_reads = []    # sr per inter-window edge that carries sr
    # Directed inter-window edges with span, for best-edge simulation.
    # Keyed by (srcWin, dstWin) -> max span-product seen (aggregate parallel
    # edges between the same window pair to their strongest span).
    inter_edge_span = {}     # (sw, dw) -> max product
    inter_edge_crossesStrand = {}  # (sw, dw) -> True if any fw<->rc

    re_wn = re.compile(r"\bwn:i:(\d+)")
    re_ws = re.compile(r"\bws:Z:(fw|rc)")
    re_sp = re.compile(r"\bsp:i:(\d+)")
    re_sn = re.compile(r"\bsn:i:(\d+)")
    re_sr = re.compile(r"\bsr:i:(\d+)")

    csv_out = None
    if csv_path is not None:
        csv_out = open(csv_path, "w")
        csv_out.write("srcWin,srcStrand,dstWin,dstStrand,sp,sn,product,sr\n")

    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                s_lines += 1
                parts = line.split("\t", 2)
                vid = parts[1]
                m_wn = re_wn.search(line)
                m_ws = re_ws.search(line)
                if m_wn:
                    wn = int(m_wn.group(1))
                    v_window[vid] = wn
                    v_strand[vid] = m_ws.group(1) if m_ws else "?"
                    win_verts[wn] = win_verts.get(wn, 0) + 1
            elif line.startswith("L\t"):
                l_lines += 1
                f = line.split("\t")
                if len(f) < 4:
                    continue
                src, dst = f[1], f[3]
                if src not in v_window or dst not in v_window:
                    continue
                sw, ss = v_window[src], v_strand[src]
                dw, ds = v_window[dst], v_strand[dst]
                edges.append((sw, ss, dw, ds))

                m_sp = re_sp.search(line)
                m_sn = re_sn.search(line)
                m_sr = re_sr.search(line)

                if sw == dw:
                    if ss == ds:
                        intra += 1
                    else:
                        inter_same_win_diff_strand += 1
                else:
                    inter_diff_win += 1
                    succ.setdefault(sw, {})
                    succ[sw][dw] = succ[sw].get(dw, 0) + 1
                    pred.setdefault(dw, {})
                    pred[dw][sw] = pred[dw].get(sw, 0) + 1

                    if m_sp and m_sn:
                        sp = int(m_sp.group(1))
                        sn = int(m_sn.group(1))
                        prod = sp * sn
                        span_products.append(prod)
                        sr = int(m_sr.group(1)) if m_sr else None
                        if sr is not None:
                            shared_reads.append(sr)
                        if csv_out:
                            csv_out.write("%s,%s,%s,%s,%d,%d,%d,%s\n" % (
                                sw, ss, dw, ds, sp, sn, prod,
                                ("" if sr is None else str(sr))))

                        key = (sw, dw)
                        if prod > inter_edge_span.get(key, -1):
                            inter_edge_span[key] = prod
                        if ss != ds:
                            inter_edge_crossesStrand[key] = True

    if csv_out:
        csv_out.close()

    # Distinct windows.
    all_wins = set(win_verts) | set(succ) | set(pred)
    n_wins = len(all_wins)

    # Tangle windows: >1 distinct predecessor OR >1 distinct successor.
    tangle_wins = multi_pred = multi_succ = multi_both = 0
    for w in all_wins:
        np_ = len(pred.get(w, {}))
        ns = len(succ.get(w, {}))
        if np_ > 1:
            multi_pred += 1
        if ns > 1:
            multi_succ += 1
        if np_ > 1 and ns > 1:
            multi_both += 1
        if np_ > 1 or ns > 1:
            tangle_wins += 1

    # Hairpin windows: window connects to its own normalized window.
    hairpin_win = set()
    for sw, ss, dw, ds in edges:
        if sw == dw and ss != ds:
            hairpin_win.add(sw)
    n_hairpin = len(hairpin_win)

    # Vertices-per-window distribution (collapse proxy).
    counts = sorted(win_verts.values())
    n_w = len(counts)
    vmin = counts[0] if n_w else 0
    vmax = counts[-1] if n_w else 0
    vsum = sum(counts)
    vmean = vsum / n_w if n_w else 0
    median = counts[n_w // 2] if n_w else 0
    p99 = counts[int(n_w * 0.99)] if n_w else 0

    if median > 0:
        collapse_thresh = 3 * median
    elif vmean > 0:
        collapse_thresh = 3 * vmean
    else:
        collapse_thresh = 1e9
    collapsed_wins = sum(1 for c in counts if c > collapse_thresh)

    print("=== Gate 1: Anchor-graph collapse/tangle stats ===")
    print("File: %s\n" % gfa_path)

    print("S lines (vertices):            %d" % s_lines)
    print("L lines (edges):               %d" % l_lines)
    print("Distinct normalized windows:   %d\n" % n_wins)

    print("-- Edge breakdown --")
    print("Intra-window (same win+strand):        %d" % intra)
    print("Strand-strand (same win, fw<->rc):     %d   <- hairpin/collapse signature"
          % inter_same_win_diff_strand)
    print("Inter-window (different window):        %d\n" % inter_diff_win)

    print("-- Window vertex counts (collapse proxy) --")
    print("min / median / mean / p99 / max:  %d / %d / %.1f / %d / %d"
          % (vmin, median, vmean, p99, vmax))
    print("Windows > 3x median verts:        %d   <- candidate collapsed windows\n"
          % collapsed_wins)

    print("-- Tangle windows --")
    print(">1 predecessor:                   %d" % multi_pred)
    print(">1 successor:                     %d" % multi_succ)
    print(">1 pred AND >1 succ:              %d" % multi_both)
    print("Tangle (>1 pred OR >1 succ):      %d  (%.1f%% of windows)"
          % (tangle_wins, (100.0 * tangle_wins / n_wins if n_wins else 0)))
    print("Hairpin windows (self fw<->rc):   %d   <- strand-strand contacts\n"
          % n_hairpin)

    # ---- Support-evidence distributions (drive N and T thresholds) ----

    print("-- Span-product (sp*sn) distribution on inter-window edges --")
    print("    (T = span-product threshold is read from here)")
    n_sp = len(span_products)
    if n_sp:
        sp = sorted(span_products)
        print("edges with sp/sn:                 %d" % n_sp)
        print("min / p10 / p25 / median:         %d / %d / %d / %d"
              % (sp[0], pct(sp, 0.10), pct(sp, 0.25), pct(sp, 0.50)))
        print("p75 / p90 / p95 / p99 / max:      %d / %d / %d / %d / %d"
              % (pct(sp, 0.75), pct(sp, 0.90), pct(sp, 0.95), pct(sp, 0.99), sp[-1]))
        print("mean:                             %.1f" % (sum(sp) / n_sp))

        hist = {}
        for v in sp:
            b = log2bucket(v)
            hist[b] = hist.get(b, 0) + 1
        low_of = {b: int(b.split("-")[0]) for b in hist}
        print("  span-product histogram (log2 buckets):")
        for b in sorted(hist, key=lambda x: low_of[x]):
            c = hist[b]
            bar = "#" * int(60 * c / n_sp + 0.5)
            print("    %-12s %8d  %s" % (b, c, bar))
    else:
        print("  (no inter-window edges carried sp/sn tags)")
    print("")

    print("-- Shared-read count (sr) distribution on inter-window edges --")
    print("    (N = minimum corroborating-read count is read from here)")
    n_sr = len(shared_reads)
    if n_sr:
        sr_hist = {}
        for v in shared_reads:
            k = "4+" if v >= 4 else str(v)
            sr_hist[k] = sr_hist.get(k, 0) + 1
        for k in ("1", "2", "3", "4+"):
            if k not in sr_hist:
                continue
            c = sr_hist[k]
            bar = "#" * int(60 * c / n_sr + 0.5)
            print("    sr=%-3s %8d  (%.1f%%)  %s"
                  % (k, c, 100.0 * c / n_sr, bar))
        print("  sr=1 edges are single-read joins (chimera-prone, count<2).")
    else:
        print("  (no inter-window edges carried sr tags)")

    # ---- Best-edge ("keep only biggest-span in/out") simulation ----
    # Simulates the destructive proposal: for each window keep only its single
    # biggest-span incoming and biggest-span outgoing inter-window edge, drop
    # the rest. Reports what that would DELETE — especially cross-strand
    # (IR fold-back) edges and edges at branch windows — so we can decide
    # whether the destructive form is safe vs the additive reciprocal-best form.
    print("\n-- Best-edge simulation (keep only biggest-span in/out per window) --")
    print("    (evaluates the destructive 'keep biggest, drop rest' proposal)")
    if inter_edge_span:
        # Best outgoing per source window, best incoming per dest window.
        best_out = {}  # sw -> (dw, span)
        best_in = {}   # dw -> (sw, span)
        for (sw, dw), span in inter_edge_span.items():
            if sw not in best_out or span > best_out[sw][1]:
                best_out[sw] = (dw, span)
            if dw not in best_in or span > best_in[dw][1]:
                best_in[dw] = (sw, span)

        total_edges = len(inter_edge_span)
        kept = set()
        for sw, (dw, _) in best_out.items():
            kept.add((sw, dw))
        for dw, (sw, _) in best_in.items():
            kept.add((sw, dw))

        dropped = [e for e in inter_edge_span if e not in kept]
        n_dropped = len(dropped)

        # How many dropped edges are cross-strand (IR fold-back risk)?
        cross_total = sum(1 for e in inter_edge_span if inter_edge_crossesStrand.get(e))
        cross_dropped = sum(1 for e in dropped if inter_edge_crossesStrand.get(e))

        # Reciprocal-best edges (the additive, safe skeleton).
        reciprocal = 0
        for sw, (dw, _) in best_out.items():
            if best_in.get(dw, (None,))[0] == sw:
                reciprocal += 1

        print("distinct directed window-pair edges:  %d" % total_edges)
        print("kept by best-in/best-out:             %d  (%.1f%%)"
              % (len(kept), 100.0 * len(kept) / total_edges))
        print("DROPPED:                              %d  (%.1f%%)"
              % (n_dropped, 100.0 * n_dropped / total_edges))
        print("  of which cross-strand (IR risk):    %d   <- destructive form would delete these"
              % cross_dropped)
        print("cross-strand edges total:             %d" % cross_total)
        if cross_total:
            print("  cross-strand edges DROPPED:         %.1f%%   <- IR loss if destructive"
                  % (100.0 * cross_dropped / cross_total))
        print("reciprocal-best edges (safe skeleton): %d  (%.1f%% of edges)"
              % (reciprocal, 100.0 * reciprocal / total_edges))
        print("  -> additive form tags these and keeps everything else (no IR loss).")
    else:
        print("  (no inter-window edges with span tags to simulate)")

    if csv_path is not None:
        print("\nPer-edge CSV written to: %s" % csv_path)


if __name__ == "__main__":
    main()
