import json, glob, math, sys, schemes
from collections import Counter
# LPB-family check: which of {sat, non-sat with the scheme's a1} reproduces KC's printed estimate
def est(Lx, Sx, Vx, a1t, t, sat):
    P = [s / l for s, l in zip(Sx, Lx)]; Q = [v / l for v, l in zip(Vx, Lx)]
    A = [0.0] * 3; B = [0.0] * 3
    for i in range(3):
        a1, a2 = 1 - 2 * P[i] - Q[i], 1 - 2 * Q[i]
        if i == t:
            if sat: continue
            a1 = a1t
        A[i] = 0.5 * math.log(1 / a1) - 0.25 * math.log(1 / a2); B[i] = 0.5 * math.log(1 / a2)
    L0, L2, L4 = Lx
    Ks = (L2 * A[1] + L4 * A[2]) / (L2 + L4) + B[2]
    Ka = A[0] + (L0 * B[0] + L2 * B[1]) / (L0 + L2)
    return Ka, Ks
res = Counter(); bad = []
pat = sys.argv[1]; methods = sys.argv[2].split(",")
for f in sorted(glob.glob(pat)):
    t = int(f.split("_")[3])
    for case in json.load(open(f)):
        for m in methods:
            r = case["kc"][m]
            val = float(r["Ka"] if t == 0 else r["Ks"]) if (r["Ka"] if t == 0 else r["Ks"]) != "NA" else 0.0
            L, S, V = schemes.accumulate(case["code"], case["s1"], case["s2"], m, {"name": "pair_count_div"})
            a1 = 1 - 2 * (S[t] / L[t]) - V[t] / L[t]
            Ka, Ks = est(L, S, V, a1, t, a1 <= 0)
            pv = Ka if t == 0 else Ks
            ok = abs(pv - val) <= 0.02 * max(1, abs(val))
            res[(m, a1 <= 0, ok)] += 1
            if not ok and len(bad) < 10: bad.append((case["name"], m, val, pv, a1))
for k in sorted(res): print(k, res[k])
for b in bad: print(b)
