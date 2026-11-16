import kc, model, itertools, random
from fractions import Fraction as F
random.seed(5)
TSU = {frozenset("AG"), frozenset("CT"), frozenset("CU")}
def tr(c): return c.replace("U", "T")
def cls(aa, c, p): return model.cls(aa, tr(c), p)
def pair(code, x, y, method):
    aa = model.T[code]; modified = method in ("MLWL", "MLPB")
    if aa[model.IDX[tr(x)]] == "*" or aa[model.IDX[tr(y)]] == "*": return [F(0)]*3, [F(0)]*3, [F(0)]*3
    L = [F(a + b, 2) for a, b in zip(model.Lc(aa, tr(x)), model.Lc(aa, tr(y)))]
    S = [F(0)]*3; V = [F(0)]*3
    diff = [p for p in range(3) if x[p] != y[p]]
    paths = []
    for order in itertools.permutations(diff):
        u = x; ok = True; st = []
        for p in order:
            v = u[:p] + y[p] + u[p+1:]
            if aa[model.IDX[tr(v)]] == "*": ok = False; break
            st.append((u, v, p)); u = v
        if ok: paths.append(st)
    for st in paths:
        for u, v, p in st:
            ts = frozenset((u[p], v[p])) in TSU
            cu, cv = cls(aa, u, p), cls(aa, v, p)
            if modified and (u, v) in model.ARG | model.ILE_S: tu = tv = True
            elif modified and (u, v) in model.ILE_V: tu = tv = False
            elif (not modified) and (u, v) in model.ARG: tu, tv = True, False
            else: tu = tv = ts
            (S if tu else V)[cu] += F(1, 2*len(paths)); (S if tv else V)[cv] += F(1, 2*len(paths))
    return L, S, V
BG1 = "GCTGGCCCAGTGTTTAAAGATCAT" * 5; BG2 = "GCC" + BG1[3:]
al = "TCAGU"
cods = ["".join(c) for c in itertools.product(al, repeat=3) if "U" in "".join(c)]
pairs = [("bg", BG1, BG2)]
for x in cods:
    for y in random.sample(kc.CODONS + cods, 40):
        pairs.append((f"{x}_{y}", x + BG1, y + BG2)); pairs.append((f"{y}_{x}", y + BG1, x + BG2))
pairs = list({n: (n, a, b) for n, a, b in pairs}.values())
bad = 0; n = 0
for m in ("LWL", "MLWL"):
    r = kc.run(pairs, methods=(m,))
    bg = kc.counts(r[("bg", m)])
    for nm, a, b in pairs[1:]:
        c = kc.counts(r[(nm, m)]); d = [[u - v for u, v in zip(p, q)] for p, q in zip(c, bg)]
        x, y = nm.split("_")
        pr = pair(1, x, y, m); n += 1
        if any(abs(float(q) - g) > 2e-4 for P, G in zip(pr, d) for q, g in zip(P, G)):
            bad += 1
            if bad < 8: print(m, x, y, "KC", [[round(v, 4) for v in z] for z in d], "model", [[float(v) for v in z] for z in pr])
print("checked", n, "bad", bad)
