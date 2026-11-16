import json, glob, math, sys, schemes
from collections import Counter
def implied(case, t):
    """from KC's printed estimate: 'sat' or ln(1/a1) of the target class."""
    m = "LWL"; r = case["kc"][m]
    L = [float(v) for v in r["Fold-Sites(0:2:4)"].split(":")]
    fold = lambda s: [0.0] * 3 if s == "NA" else [float(v) for v in s.split(":")]
    Si, Vi = fold(r["Fold-S-Substitutions(0:2:4)"]), fold(r["Fold-N-Substitutions(0:2:4)"])
    # exact counts from the model (rational) for the other terms
    Lx, Sx, Vx = schemes.accumulate(case["code"], case["s1"], case["s2"], m, {"name": "pair_count_div"})
    P = [s / l for s, l in zip(Sx, Lx)]; Q = [v / l for v, l in zip(Vx, Lx)]
    A = [0.0] * 3; B = [0.0] * 3
    for i in range(3):
        if i == t: continue
        a1, a2 = 1 - 2 * P[i] - Q[i], 1 - 2 * Q[i]
        A[i] = max(0.0, 0.5 * math.log(1 / a1) - 0.25 * math.log(1 / a2)); B[i] = 0.5 * math.log(1 / a2)
    bt = 1 - 2 * Q[t]; Bt = 0.5 * math.log(1 / bt)
    L0, L2, L4 = Lx
    if t == 0:
        val = float(r["Ka"]) if r["Ka"] != "NA" else 0.0
        # Ka = 3 (L2 B2 + L0 K0) / (2 L2 + 3 L0), K0 = A0 + B0; sat -> K0 = 0
        sat = 3 * (L2 * B[1]) / (2 * L2 + 3 * L0)
        K0 = val * (2 * L2 + 3 * L0) / 3 / L0 - L2 * B[1] / L0
        A0 = K0 - Bt
        x = (A0 + 0.25 * math.log(1 / bt)) / 0.5
    elif t == 1:
        val = float(r["Ks"]) if r["Ks"] != "NA" else 0.0
        sat = 3 * (L4 * (A[2] + B[2])) / (L2 + 3 * L4)
        A2 = (val * (L2 + 3 * L4) / 3 - L4 * (A[2] + B[2])) / L2
        x = (A2 + 0.25 * math.log(1 / bt)) / 0.5
    else:
        val = float(r["Ks"]) if r["Ks"] != "NA" else 0.0
        sat = 3 * (L2 * A[1]) / (L2 + 3 * L4)
        K4 = (val * (L2 + 3 * L4) / 3 - L2 * A[1]) / L4
        A4 = K4 - Bt
        x = (A4 + 0.25 * math.log(1 / bt)) / 0.5
    if abs(val - sat) <= 1e-5 * max(1, abs(sat)):
        return "sat", val
    return x, val

res = Counter(); rows = []
for f in sorted(glob.glob("zero_LWL_*.json")):
    t = int(f.split("_")[3])
    for case in json.load(open(f)):
        imp, val = implied(case, t)
        row = {"name": case["name"], "kc": imp}
        for sc in schemes.SCHEMES:
            L, S, V = schemes.accumulate(case["code"], case["s1"], case["s2"], "LWL", sc)
            a1 = 1 - 2 * (S[t] / L[t]) - V[t] / L[t]
            pred = "sat" if a1 <= 0 else math.log(1 / a1)
            ok = (pred == "sat" and imp == "sat") or (pred != "sat" and imp != "sat" and abs(pred - imp) < 1e-3 * abs(pred))
            res[(sc["name"], ok)] += 1
            row[sc["name"]] = pred
        rows.append(row)
for k in sorted(res): print(k, res[k])
print(Counter(r["kc"] == "sat" for r in rows))
for r in rows[:8]: print(r)
print("--- misses of pair_count_div")
for r in rows:
    p = r["pair_count_div"]; k = r["kc"]
    ok = (p == "sat" and k == "sat") or (p != "sat" and k != "sat" and abs(p - k) < 1e-3 * abs(p))
    if not ok: print(r["name"], k, p, r["pair_count_mulw"], r["pair_acc_half"])
