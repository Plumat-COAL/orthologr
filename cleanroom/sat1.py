import kc, math
# 4-fold class saturated: GCT->GCC (S4) many, GCT/GCT few; plus 2-fold and 0-fold material
a = "GCT"*6 + "TTTAAAGATCAT"*5 + "ATGTGG"*3
b = "GCC"*5 + "GCA" + "TTCAAGGACCAC"*1 + "TTTAAAGATCAT"*4 + "ATGTGG"*3
b = b
r = kc.run([("s", a, b)])
for (n, m), row in r.items():
    print(m, row["Ka"], row["Ks"], row["Fold-Sites(0:2:4)"], row["Fold-S-Substitutions(0:2:4)"], row["Fold-N-Substitutions(0:2:4)"])
