import sys, json, kc
from fractions import Fraction
BG1 = "GCTGGCCCAGTGTTTAAAGATCAT" * 5
BG2 = "GCC" + BG1[3:]
out = {}
for code in kc.CODES:
    pairs = [("bg", BG1, BG2)] + [(f"{x}_{y}", x + BG1, y + BG2) for x in kc.CODONS for y in kc.CODONS]
    r = kc.run(pairs, code=code)
    d = {}
    for (n, m), row in r.items():
        d[f"{n}|{m}"] = kc.counts(row)
    out[code] = d
    print(code, len(d), file=sys.stderr)
json.dump(out, open("pairs.json", "w"))
