import json, kc, model, sys
from collections import Counter
P = json.load(open("pairs.json"))
def sub(a, b): return [[round(x - y, 4) for x, y in zip(u, v)] for u, v in zip(a, b)]
tot = Counter(); ex = []
codes = [int(c) for c in sys.argv[1:]] or kc.CODES
for code in codes:
    d = P[str(code)]
    for m in ["LWL", "MLWL", "LPB", "MLPB"]:
        bg = d[f"bg|{m}"]
        for x in kc.CODONS:
            for y in kc.CODONS:
                if x == y: continue
                r = sub(d[f"{x}_{y}|{m}"], bg)
                pr = [[round(float(v), 4) for v in vec] for vec in model.pair(code, x, y, m)]
                k = sum(a != b for a, b in zip(x, y))
                ok = all(abs(a - b) < 2e-4 for u, v in zip(r, pr) for a, b in zip(u, v))
                tot[(k, ok)] += 1
                if not ok and len(ex) < 30: ex.append((code, m, x, y, "KC", r, "model", pr))
print(tot)
for e in ex: print(*e)
