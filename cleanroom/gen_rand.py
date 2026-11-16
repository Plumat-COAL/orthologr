"""Random/simulated pairs for all codes; KC output -> TSV for the R comparison."""
import random, kc, sys, csv
random.seed(int(sys.argv[1])); n = int(sys.argv[2]); out = sys.argv[3]
alpha = "TCAG"
rows = []
for code in kc.CODES:
    pairs = []
    for i in range(n):
        L = random.randint(1, 400)
        s1 = [random.choice(kc.CODONS) for _ in range(L)]
        div = random.choice([0.02, 0.1, 0.3, 0.6, 1.0])
        s2 = []
        for c in s1:
            c = list(c)
            for p in range(3):
                if random.random() < div * random.random():
                    c[p] = random.choice(alpha)
            s2.append("".join(c))
        a, b = "".join(s1), "".join(s2)
        # sprinkle odd characters / lowercase / U in some pairs
        r = random.random()
        if r < 0.15:
            a = "".join(ch if random.random() > 0.02 else random.choice("-NnRU") for ch in a)
        elif r < 0.25:
            b = b.lower()
        elif r < 0.3:
            b = b.replace("T", "U")
        pairs.append((f"r{code}_{i}", a, b))
    res = kc.run(pairs, code=code)
    for nm, a, b in pairs:
        for m in ("LWL", "MLWL", "LPB", "MLPB"):
            row = res.get((nm, m))
            if row is None:
                rows.append([nm, code, m, a, b, "MISSING", "", ""]); continue
            rows.append([nm, code, m, a, b, row["Fold-Sites(0:2:4)"], row["Fold-S-Substitutions(0:2:4)"], row["Fold-N-Substitutions(0:2:4)"]])
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["name", "code", "method", "s1", "s2", "L", "S", "V"]); w.writerows(rows)
print(len(rows), sum(r[5] == "MISSING" for r in rows))
