#!/usr/bin/env python3
"""Score sites against the chr12 truth track, partitioned on whether the TRUTH
site sits in a homopolymer.

The truth track is itself a minimap2 mat-vs-pat alignment, so inside a
homopolymer the placement of a difference is ambiguous -- the same change spells
as a substitution or as a shifted indel. Such a call can neither confirm nor
refute a detector, so scoring a homopolymer GATE against it measures agreement
between two guesses. This reports both partitions so that trap is visible
rather than averaged away.
"""
"""Precision/recall for a sites.tsv against the chr12 mat-vs-pat het SNV truth.

Projects each site (owner read, read position) into hg002v1.1 coordinates by
walking an independent minimap2 alignment of that read, then matches in the
frame the read aligned to. Absolute numbers depend on the coverage denominator;
what this is for is COMPARING two site sets on identical footing.
"""
import re, sys
from collections import Counter, defaultdict
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# Reference, to ask of each TRUTH site whether it sits in a homopolymer.
ref = {}
_name = None
for _line in open('HG002_chr12:11000000-17000000.fasta'):
    if _line.startswith('>'):
        _name = _line[1:].split()[0].split(':')[0]; ref[_name] = []
    else:
        ref[_name].append(_line.strip())
ref = {k: ''.join(v).upper() for k, v in ref.items()}
REFOFF = 11000000

def in_homopolymer(contig, pos1):
    "dinara's own unit-1 test: 3 identical bases immediately before or after."
    s_ = ref.get(contig)
    if s_ is None: return False
    i = pos1 - 1 - REFOFF
    if not (3 <= i < len(s_) - 3): return False
    after, before = s_[i+1:i+4], s_[i-3:i]
    return len(set(after)) == 1 or len(set(before)) == 1

truth = {}          # (contig, pos1) -> canonical id (maternal pos)
truth_hp = {}       # canonical id -> bool
canon = set()
with open('truth_chr12_11_17_snv.tsv') as fh:
    fh.readline()
    for line in fh:
        mc, mp, pc, pp, ma, pa = line.rstrip('\n').split('\t')
        if len(ma) != 1 or len(pa) != 1:
            continue
        truth[(mc, int(mp))] = int(mp)
        truth[(pc, int(pp))] = int(mp)
        canon.add(int(mp))
        truth_hp[int(mp)] = (in_homopolymer(mc, int(mp)) or
                             in_homopolymer(pc, int(pp)))

# The truth track covers only part of the requested window, so a site outside its
# span cannot be scored either way -- counting those as false positives just
# measures where the track stops. Restrict both numerator and denominator to it.
truth_span = {}
for (c, p_) in truth:
    lo, hi = truth_span.get(c, (p_, p_))
    truth_span[c] = (min(lo, p_), max(hi, p_))
print("truth span: " + str(truth_span))

best = {}
for line in open('reads_vs_chr12.paf'):
    f = line.rstrip('\n').split('\t')
    cigar = next((t[5:] for t in f[12:] if t.startswith('cg:Z:')), None)
    if cigar is None:
        continue
    # Contigs are named "<hap>:11000000-17000000" and positions are relative to
    # that window; the truth track is absolute on "<hap>". Normalise both here.
    tname = f[5]
    offset = 0
    if ':' in tname:
        tname, span = tname.split(':', 1)
        offset = int(span.split('-')[0])
    q = dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
             tname=tname, ts=int(f[7]) + offset, te=int(f[8]) + offset,
             nmatch=int(f[9]), cigar=cigar)
    if f[0] not in best or q['nmatch'] > best[f[0]]['nmatch']:
        best[f[0]] = q

def projector(a):
    m = {}
    q = a['qs'] if a['strand'] == '+' else (a['qlen'] - a['qe'])
    t = a['ts']
    for n, op in CIGAR_RE.findall(a['cigar']):
        n = int(n)
        if op in 'M=X':
            for i in range(n):
                orig = (q + i) if a['strand'] == '+' else (a['qlen'] - 1 - (q + i))
                m[orig] = t + i + 1
            q += n; t += n
        elif op in 'IS':
            q += n
        elif op in 'DN':
            t += n
    return m

def score(path, tol=2):
    sites = []
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        sites.append(dict(read=f[0], pos=int(f[1])))
    by_read = defaultdict(list)
    for s in sites:
        by_read[s['read']].append(s)
    projected = []
    for read, ss in by_read.items():
        a = best.get(read)
        if a is None:
            continue
        m = projector(a)
        for s in ss:
            tp = m.get(s['pos'])
            if tp is None:
                continue
            lo_hi = truth_span.get(a["tname"])
            if lo_hi is None or not (lo_hi[0] <= tp <= lo_hi[1]):
                continue
            s['tname'] = a['tname']; s['tpos'] = tp
            projected.append(s)
    def match(s):
        for d in range(-tol, tol + 1):
            c = truth.get((s['tname'], s['tpos'] + d))
            if c is not None:
                return c
        return None
    hits, found = 0, set()
    for s in projected:
        c = match(s)
        if c is not None:
            hits += 1; found.add(c)
    found_nhp = sum(1 for c in found if not truth_hp.get(c, False))
    return len(sites), len(projected), hits, len(found), found_nhp

print(f"truth: {len(canon)} distinct het SNVs")
nhp_total = sum(1 for c in canon if not truth_hp.get(c, False))
hp_total = len(canon) - nhp_total
print(f"truth partition: {nhp_total} non-homopolymer, {hp_total} homopolymer (unverifiable)")
print(f"{'set':20s} {'sites':>7s} {'FP':>5s} {'precision':>10s} {'found(all)':>11s} {'found(non-HP)':>14s} {'recall(non-HP)':>15s}")
for path in sys.argv[1:]:
    n, proj, hits, found, found_nhp = score(path)
    prec = 100.0 * hits / proj if proj else 0.0
    print(f"{path:20s} {n:7d} {proj-hits:5d} {prec:9.1f}% {found:11d} {found_nhp:14d} "
          f"{100.0*found_nhp/nhp_total:14.1f}%")
