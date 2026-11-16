import itertools, kc, ncbi
from fractions import Fraction as F
T = ncbi.tables()
IDX = {x: i for i, x in enumerate(kc.CODONS)}
TS = {frozenset("AG"), frozenset("CT")}
ARG = {("CGA", "AGA"), ("AGA", "CGA"), ("CGG", "AGG"), ("AGG", "CGG")}
ILE_S = {("ATA", "ATT"), ("ATT", "ATA"), ("ATA", "ATC"), ("ATC", "ATA")}
ILE_V = {("ATA", "ATG"), ("ATG", "ATA")}

def cls(aa, c, p):
    a = aa[IDX[c]]
    syn = sum(aa[IDX[c[:p] + n + c[p + 1:]]] == a for n in "TCAG" if n != c[p])
    return 0 if syn == 0 else (2 if syn == 3 else 1)

def Lc(aa, c):
    L = [0, 0, 0]
    if aa[IDX[c]] == "*":
        return L
    for p in range(3):
        L[cls(aa, c, p)] += 1
    return L

def step(aa, u, v, p, modified):
    """(S, V) vectors (Fractions) for the single step u -> v at position p."""
    S = [F(0)] * 3; V = [F(0)] * 3
    ts = frozenset((u[p], v[p])) in TS
    cu, cv = cls(aa, u, p), cls(aa, v, p)
    if modified and (u, v) in ARG | ILE_S:
        tu = tv = True
    elif modified and (u, v) in ILE_V:
        tu = tv = False
    elif (not modified) and (u, v) in ARG:
        tu, tv = True, False
    else:
        tu = tv = ts
    (S if tu else V)[cu] += F(1, 2)
    (S if tv else V)[cv] += F(1, 2)
    return S, V

def pair(code, x, y, method):
    aa = T[code]; modified = method in ("MLWL", "MLPB")
    S = [F(0)] * 3; V = [F(0)] * 3
    if aa[IDX[x]] == "*" or aa[IDX[y]] == "*":
        return [F(0)] * 3, S, V
    L = [F(a + b, 2) for a, b in zip(Lc(aa, x), Lc(aa, y))]
    diff = [p for p in range(3) if x[p] != y[p]]
    paths = []
    for order in itertools.permutations(diff):
        u = x; ok = True; steps = []
        for p in order:
            v = u[:p] + y[p] + u[p + 1:]
            if aa[IDX[v]] == "*":
                ok = False; break
            steps.append((u, v, p)); u = v
        if ok:
            paths.append(steps)
    if not paths:
        return L, S, V
    w = F(1, len(paths))
    for steps in paths:
        for u, v, p in steps:
            s, t = step(aa, u, v, p, modified)
            S = [a + w * b for a, b in zip(S, s)]
            V = [a + w * b for a, b in zip(V, t)]
    return L, S, V
