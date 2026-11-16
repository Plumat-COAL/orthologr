import itertools, model, kc
from fractions import Fraction as F

def steps_halves(code, x, y, method):
    """list of paths; each path = list of (S/V, class) halves, in step order (u-half then v-half)."""
    aa = model.T[code]; modified = method in ("MLWL", "MLPB")
    if aa[model.IDX[x]] == "*" or aa[model.IDX[y]] == "*":
        return None
    diff = [p for p in range(3) if x[p] != y[p]]
    paths = []
    for order in itertools.permutations(diff):
        u = x; ok = True; halves = []
        for p in order:
            v = u[:p] + y[p] + u[p + 1:]
            if aa[model.IDX[v]] == "*":
                ok = False; break
            s, t = model.step(aa, u, v, p, modified)
            cu, cv = model.cls(aa, u, p), model.cls(aa, v, p)
            # recover the two halves' kinds
            ku = "S" if s[cu] >= F(1, 2) or (cu == cv and s[cu] == 1) else "V"
            # recompute explicitly
            ts = frozenset((u[p], v[p])) in model.TS
            if modified and (u, v) in model.ARG | model.ILE_S: tu = tv = True
            elif modified and (u, v) in model.ILE_V: tu = tv = False
            elif (not modified) and (u, v) in model.ARG: tu, tv = True, False
            else: tu = tv = ts
            halves.append(("S" if tu else "V", cu)); halves.append(("S" if tv else "V", cv))
            u = v
        if ok:
            paths.append(halves)
    return paths

def Lpair(code, x, y):
    aa = model.T[code]
    if aa[model.IDX[x]] == "*" or aa[model.IDX[y]] == "*":
        return [0.0, 0.0, 0.0]
    return [ (a + b) / 2 for a, b in zip(model.Lc(aa, x), model.Lc(aa, y))]

def accumulate(code, s1, s2, method, scheme):
    L = [0.0] * 3; S = [0.0] * 3; V = [0.0] * 3
    for i in range(0, len(s1), 3):
        x, y = s1[i:i + 3], s2[i:i + 3]
        l = Lpair(code, x, y)
        if scheme.get("Lhalf"):
            aa = model.T[code]
            if aa[model.IDX[x]] != "*" and aa[model.IDX[y]] != "*":
                for c in (model.Lc(aa, x), model.Lc(aa, y)):
                    for k in range(3): L[k] += c[k] * 0.5
        else:
            for k in range(3): L[k] += l[k]
        if x == y: continue
        paths = steps_halves(code, x, y, method)
        if not paths: continue
        n = len(paths); sc = scheme["name"]
        if sc == "global_half":          # tot += 0.5/n for each half, as they come
            w = 0.5 / n
            for h in paths:
                for kind, c in h: (S if kind == "S" else V)[c] += w
        elif sc == "global_w_times_half":  # tot += (1/n) * 0.5
            w = 1.0 / n
            for h in paths:
                for kind, c in h: (S if kind == "S" else V)[c] += w * 0.5
        elif sc == "pair_acc_half":       # pair accumulators of 0.5/n, then add
            ps = [0.0] * 3; pv = [0.0] * 3; w = 0.5 / n
            for h in paths:
                for kind, c in h: (ps if kind == "S" else pv)[c] += w
            for k in range(3): S[k] += ps[k]; V[k] += pv[k]
        elif sc == "pair_count_div":      # exact half counts / n
            ps = [0.0] * 3; pv = [0.0] * 3
            for h in paths:
                for kind, c in h: (ps if kind == "S" else pv)[c] += 0.5
            for k in range(3): S[k] += ps[k] / n; V[k] += pv[k] / n
        elif sc == "pair_count_mulw":
            ps = [0.0] * 3; pv = [0.0] * 3
            for h in paths:
                for kind, c in h: (ps if kind == "S" else pv)[c] += 0.5
            w = 1.0 / n
            for k in range(3): S[k] += ps[k] * w; V[k] += pv[k] * w
        elif sc == "path_acc":            # per-path exact sums, pair += path/n
            ps = [0.0] * 3; pv = [0.0] * 3
            for h in paths:
                qs = [0.0] * 3; qv = [0.0] * 3
                for kind, c in h: (qs if kind == "S" else qv)[c] += 0.5
                for k in range(3): ps[k] += qs[k] / n; pv[k] += qv[k] / n
            for k in range(3): S[k] += ps[k]; V[k] += pv[k]
        elif sc == "path_global":         # per-path exact sums, tot += path/n
            for h in paths:
                qs = [0.0] * 3; qv = [0.0] * 3
                for kind, c in h: (qs if kind == "S" else qv)[c] += 0.5
                for k in range(3): S[k] += qs[k] / n; V[k] += qv[k] / n
        elif sc == "step_global":         # per step: both halves summed (exact) then * 1/n ... tot += stepsum/n
            for h in paths:
                for j in range(0, len(h), 2):
                    qs = [0.0] * 3; qv = [0.0] * 3
                    for kind, c in h[j:j + 2]: (qs if kind == "S" else qv)[c] += 0.5
                    for k in range(3):
                        if qs[k]: S[k] += qs[k] / n
                        if qv[k]: V[k] += qv[k] / n
        else:
            raise ValueError(sc)
    return L, S, V

SCHEMES = [{"name": n} for n in ("global_half", "global_w_times_half", "pair_acc_half", "pair_count_div", "pair_count_mulw", "path_acc", "path_global", "step_global")]
