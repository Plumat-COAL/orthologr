import kc, model
BG1 = "GCTGGCCCAGTGTTTAAAGATCAT" * 5
BG2 = "GCC" + BG1[3:]
odd = ["---", "GC-", "-CT", "NNN", "GCN", "NCT", "gct", "gcc", "GCU", "UUU", "GCR", "ttt", "XXX", "G-T", "?CT", "TAA", "taa"]
targets = ["GCT", "GCC", "TTT", "TTC", "ATG", "TAA", "---", "GCN", "gct", "TCT"]
pairs = [("bg", BG1, BG2)]
for x in odd:
    for y in targets:
        pairs.append((f"{x}|{y}", x + BG1, y + BG2))
        pairs.append((f"{y}|{x}", y + BG1, x + BG2))
r = kc.run(pairs, methods=("LWL",))
bg = kc.counts(r[("bg", "LWL")])
seen = set()
for n, a, b in pairs[1:]:
    if n in seen: continue
    seen.add(n)
    row = r.get((n, "LWL"))
    if row is None: print(n, "NO ROW"); continue
    c = kc.counts(row)
    d = [[round(u - v, 4) for u, v in zip(p, q)] for p, q in zip(c, bg)]
    x, y = n.split("|")
    up = lambda s: s.upper().replace("U", "T")
    try:
        pr = [[float(v) for v in vec] for vec in model.pair(1, up(x), up(y), "LWL")]
    except Exception:
        pr = "n/a"
    print(f"{n:10s} len={row['Length']:4s} KC={d} model(upper,U->T)={pr}")
