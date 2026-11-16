import random, json, kc, model, sys
from fractions import Fraction as F
seed, target, ncase = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
method = sys.argv[4] if len(sys.argv) > 4 else "LWL"
code = int(sys.argv[5]) if len(sys.argv) > 5 else 1
random.seed(seed); aa = model.T[code]
sense = [c for c in kc.CODONS if aa[model.IDX[c]] != "*"]
C = {(x, y): model.pair(code, x, y, method) for x in sense for y in sense}
def dlt(c, t): L, S, V = c; return L[t] - 2 * S[t] - V[t]
frac = [k for k, c in C.items() if (c[1][target].denominator > 1 or c[2][target].denominator > 1)]
integral = [k for k, c in C.items() if all(v.denominator == 1 for v in c[1] + c[2]) and dlt(c, target).denominator == 1]
plus = [k for k in integral if dlt(C[k], target) == 1]
minus = [k for k in integral if dlt(C[k], target) == -1]
if not minus:
    half = [k for k, c in C.items() if dlt(c, target) == F(-1, 2) and all(v.denominator == 1 for v in c[1] + c[2])]
    minus = [None]
filler = [k for k, c in C.items() if k[0] == k[1]]
assert plus and minus, (plus, minus)
cases = []
from collections import Counter
REJ = Counter(); it = 0
while len(cases) < ncase:
    it += 1
    if it % 2000 == 0: print(it, len(cases), REJ, file=sys.stderr, flush=True)
    if it > 20000: break
    sel = [random.choice(frac) for _ in range(random.randint(2, 30))]
    D = sum(dlt(C[k], target) for k in sel)
    n = 0
    while (D - (D.numerator // D.denominator)) != 0 and n < 40:
        k = random.choice(frac); sel.append(k); D += dlt(C[k], target); n += 1
    if (D - (D.numerator // D.denominator)) != 0:
        REJ['frac'] += 1; continue
    fill = [random.choice(filler) for _ in range(random.randint(20, 120))]
    D = int(D + sum(dlt(C[k], target) for k in fill))
    if D > 0 and minus == [None]:
        extra = [random.choice(half) for _ in range(2 * D)]
    else:
        extra = [random.choice(minus) for _ in range(D)] if D > 0 else [random.choice(plus) for _ in range(-D)]
    allp = sel + extra + fill
    random.shuffle(allp)
    L = [sum(C[k][0][i] for k in allp) for i in range(3)]
    S = [sum(C[k][1][i] for k in allp) for i in range(3)]
    V = [sum(C[k][2][i] for k in allp) for i in range(3)]
    if L[target] - 2 * S[target] - V[target] != 0 or 1 - 2 * V[target] / L[target] <= 0:
        REJ['t'] += 1; continue
    if not all(L[i] > 0 and 1 - 2 * S[i] / L[i] - V[i] / L[i] > F(1, 50) and 1 - 2 * V[i] / L[i] > F(1, 50) for i in range(3) if i != target):
        REJ['o'] += 1; continue
    cases.append(("".join(x for x, _ in allp), "".join(y for _, y in allp)))
pairs = [(f"z{target}_{seed}_{i}", a, b) for i, (a, b) in enumerate(cases)]
r = kc.run(pairs, methods=("LWL", "LPB", "MLWL", "MLPB"), code=code)
keys = ("Ka", "Ks", "Fold-Sites(0:2:4)", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)")
out = [{"name": n, "code": code, "s1": a, "s2": b, "kc": {m: {k: r[(n, m)][k] for k in keys} for m in ("LWL", "LPB", "MLWL", "MLPB")}} for n, a, b in pairs]
json.dump(out, open(f"zero_{method}_{code}_{target}_{seed}.json", "w"))
print(len(out))
