#!/usr/bin/env python3
"""Dump and characterise the false positives of a sites.tsv.

A site is a FP when no truth het SNV sits within +/-tol of its projected
position. That is not the same as "wrong": the truth track is itself a
minimap2 mat-vs-pat alignment, so inside a homopolymer or short tandem repeat
the placement of a difference is ambiguous and a call there can neither be
confirmed nor refuted. This splits the FPs on that, and reports how far the
nearest truth variant actually is.
"""
import re, sys
from collections import Counter, defaultdict
CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

ref = {}
name = None
for line in open('HG002_chr12:11000000-17000000.fasta'):
    if line.startswith('>'):
        name = line[1:].split()[0].split(':')[0]
        ref[name] = []
    else:
        ref[name].append(line.strip())
ref = {k: ''.join(v) for k, v in ref.items()}
OFFSET = 11000000

truth = {}
with open('truth_chr12_11_17_snv.tsv') as fh:
    fh.readline()
    for line in fh:
        mc, mp, pc, pp, ma, pa = line.rstrip('\n').split('\t')
        if len(ma) == 1 and len(pa) == 1:
            truth[(mc, int(mp))] = int(mp)
            truth[(pc, int(pp))] = int(mp)
truth_by_contig = defaultdict(list)
for (c, p) in truth:
    truth_by_contig[c].append(p)
for c in truth_by_contig:
    truth_by_contig[c].sort()

best = {}
for line in open('reads_vs_chr12.paf'):
    f = line.rstrip('\n').split('\t')
    cigar = next((t[5:] for t in f[12:] if t.startswith('cg:Z:')), None)
    if cigar is None: continue
    tname = f[5]; off = 0
    if ':' in tname:
        tname, span = tname.split(':', 1); off = int(span.split('-')[0])
    q = dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], tname=tname,
             ts=int(f[7])+off, te=int(f[8])+off, nmatch=int(f[9]), cigar=cigar)
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
                orig = (q+i) if a['strand']=='+' else (a['qlen']-1-(q+i))
                m[orig] = t + i + 1
            q += n; t += n
        elif op in 'IS': q += n
        elif op in 'DN': t += n
    return m

def context(contig, pos1, w=6):
    s = ref.get(contig)
    if s is None: return ''
    i = pos1 - 1 - OFFSET
    return s[max(0, i-w): i+w+1].upper() if 0 <= i < len(s) else ''

def homopolymer(ctx, w=6):
    """dinara's own test: >=3 identical bases immediately before or after."""
    if len(ctx) < w+4: return False
    before, base, after = ctx[:w], ctx[w], ctx[w+1:]
    return (after[:3] and len(set(after[:3])) == 1) or \
           (before[-3:] and len(set(before[-3:])) == 1)

def str_repeat(ctx, w=6):
    for unit in (2, 3, 4):
        seg = ctx[w+1: w+1+unit*3]
        if len(seg) == unit*3 and seg[:unit]*3 == seg: return True
        seg = ctx[max(0, w-unit*3): w]
        if len(seg) == unit*3 and seg[:unit]*3 == seg: return True
    return False

path, tol = sys.argv[1], 2
sites = []
for line in open(path):
    f = line.rstrip('\n').split('\t')
    sites.append(dict(read=f[0], pos=int(f[1]), dom=f[2], alt=f[3],
                      domc=int(f[4]), altc=int(f[5])))
by_read = defaultdict(list)
for s in sites: by_read[s['read']].append(s)
proj = []
for read, ss in by_read.items():
    a = best.get(read)
    if a is None: continue
    m = projector(a)
    for s in ss:
        tp = m.get(s['pos'])
        if tp is not None:
            s['tname'] = a['tname']; s['tpos'] = tp; proj.append(s)

fps = []
for s in proj:
    if not any(truth.get((s['tname'], s['tpos']+d)) for d in range(-tol, tol+1)):
        fps.append(s)

print(f"{path}: {len(sites)} sites, {len(proj)} projected, {len(fps)} FP at +/-{tol}")
cnt = Counter()
for s in fps:
    ctx = context(s['tname'], s['tpos'])
    hp, st = homopolymer(ctx), str_repeat(ctx)
    tl = truth_by_contig[s['tname']]
    import bisect
    i = bisect.bisect_left(tl, s['tpos'])
    near = min([abs(tl[j]-s['tpos']) for j in (i-1, i) if 0 <= j < len(tl)] or [10**9])
    cls = 'homopolymer' if hp else ('STR' if st else ('near-truth<=10bp' if near <= 10 else 'isolated'))
    cnt[cls] += 1
    s['ctx'], s['near'], s['cls'] = ctx, near, cls
for k, v in cnt.most_common():
    print(f"    {k:20s} {v}")
print("    ---- detail ----")
for s in sorted(fps, key=lambda x: x['near'])[:25]:
    print(f"    {s['tname'][-8:]}:{s['tpos']} {s['dom']}/{s['alt']} "
          f"{s['domc']}/{s['altc']}  ctx={s['ctx']}  nearestTruth={s['near']}  {s['cls']}")
