"""Black-box helpers: run the KaKs_Calculator binary on constructed AXT input, parse its printed table."""
import subprocess, os, tempfile, itertools

KC = "/home/claude/tools-g/bin/KaKs_Calculator"
CODES = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23]
B = "TCAG"
CODONS = ["".join(c) for c in itertools.product(B, repeat=3)]


def run(pairs, methods=("LWL", "MLWL", "LPB", "MLPB"), code=1, workdir="/home/claude/w-olclean-scratch/tmp"):
    """pairs: list of (name, s1, s2). Returns {(name, method): row dict}."""
    os.makedirs(workdir, exist_ok=True)
    fd, inp = tempfile.mkstemp(suffix=".axt", dir=workdir)
    with os.fdopen(fd, "w") as f:
        for n, a, b in pairs:
            f.write(f"{n}\n{a}\n{b}\n\n")
    out = inp + ".kaks"
    args = [KC, "-i", inp, "-o", out, "-c", str(code)]
    for m in methods:
        args += ["-m", m]
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    res = {}
    with open(out) as f:
        hdr = None
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[0] == "Sequence":
                hdr = x
                continue
            row = dict(zip(hdr, x))
            res[(x[0], x[1])] = row
    os.remove(inp)
    os.remove(out)
    return res


def fold(s):
    if s == "NA":
        return [0.0, 0.0, 0.0]
    return [float(v) for v in s.split(":")]


def counts(row):
    return (fold(row["Fold-Sites(0:2:4)"]), fold(row["Fold-S-Substitutions(0:2:4)"]),
            fold(row["Fold-N-Substitutions(0:2:4)"]))
