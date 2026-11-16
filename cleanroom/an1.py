import json, kc, ncbi
from fractions import Fraction
P = json.load(open("pairs.json"))
T = ncbi.tables()
def sub(a, b): return [[round(x - y, 4) for x, y in zip(u, v)] for u, v in zip(a, b)]
ts = {("A","G"),("G","A"),("C","T"),("T","C")}
def classes(code, c):
    aa = T[code]; idx = {x:i for i,x in enumerate(kc.CODONS)}
    a = aa[idx[c]]
    out = []
    for p in range(3):
        syn = 0
        for n in "TCAG":
            if n == c[p]: continue
            d = c[:p]+n+c[p+1:]
            if aa[idx[d]] == a: syn += 1
        out.append(0 if syn == 0 else (2 if syn == 4-1 and False else (4 if syn == 3 else 2)))
    return out
bad = 0
for code in kc.CODES:
    d = P[str(code)]
    for m in ["LWL","MLWL","LPB","MLPB"]:
        bg = d[f"bg|{m}"]
        for c in kc.CODONS:
            r = d.get(f"{c}_{c}|{m}")
            L = sub(r, bg)[0]
            aa = T[code][kc.CODONS.index(c)]
            pred = [0,0,0]
            if aa != "*":
                for k in classes(code, c): pred[{0:0,2:1,4:2}[k]] += 1
            if L != pred:
                bad += 1
                if bad < 40: print(code, m, c, aa, "KC", L, "pred", pred, "S/V", sub(r,bg)[1:])
print("bad", bad)
