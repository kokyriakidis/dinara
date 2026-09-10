import sys
kept = set()
for line in open(sys.argv[2]):
    r, p = line.rstrip('\n').split('\t')
    kept.add((r, int(p)))
n = 0
with open(sys.argv[3], 'w') as out:
    for line in open(sys.argv[1]):
        f = line.rstrip('\n').split('\t')
        if (f[0], int(f[1])) in kept:
            out.write(line); n += 1
print(f"{sys.argv[3]}: {n} sites")
