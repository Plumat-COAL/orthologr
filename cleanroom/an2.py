import json, kc, ncbi, sys
from collections import Counter
P = json.load(open("pairs.json"))
T = ncbi.tables()
IDX = {x:i for i,x in enumerate(kc.CODONS)}
def sub(a, b): return [[round(x - y, 4) for x, y in zip(u, v)] for u, v in zip(a, b)]
TS = {frozenset("AG"), frozenset("CT")}
def cls(aa, c, p):
    a = aa[IDX[c]]; syn = sum(aa[IDX[c[:p]+n+c[p+1:]]] == a for n in "TCAG" if n != c[p])
    return 0 if syn == 0 else (2 if syn == 3 else 1)   # column index 0/1/2 for 0,2,4-fold
code = int(sys.argv[1]); m = sys.argv[2]
aa = T[code]; d = P[str(code)]; bg = d[f"bg|{m}"]
stats = Counter(); ex = []
for x in kc.CODONS:
    for y in kc.CODONS:
        diff = [p for p in range(3) if x[p] != y[p]]
        if len(diff) != 1: continue
        r = sub(d[f"{x}_{y}|{m}"], bg)
        p = diff[0]; ts = frozenset((x[p], y[p])) in TS
        cx, cy = (cls(aa, x, p) if aa[IDX[x]] != "*" else None), (cls(aa, y, p) if aa[IDX[y]] != "*" else None)
        S, V = r[1], r[2]
        key = (aa[IDX[x]]=="*", aa[IDX[y]]=="*", cx == cy)
        # hypothesis: average of the classes in x and y
        pred = [0.0]*3
        cl = [c for c in (cx, cy) if c is not None]
        if aa[IDX[x]] != "*" and aa[IDX[y]] != "*":
            for c in cl: pred[c] += 0.5
        got = S if ts else V
        other = V if ts else S
        ok = [round(g,3) for g in got] == [round(q,3) for q in pred] and all(o == 0 for o in other)
        stats[(key, ok)] += 1
        if not ok and len(ex) < 25: ex.append((x, y, aa[IDX[x]], aa[IDX[y]], "ts" if ts else "tv", cx, cy, "S", S, "V", V, "L", r[0]))
print(stats)
for e in ex: print(*e)
