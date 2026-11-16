import re
def tables(path="/home/claude/w-olclean-scratch/gc.prt"):
    txt = open(path).read()
    t = {}
    for m in re.finditer(r'id (\d+) ,\s*ncbieaa\s+"([A-Z*]{64})"', txt):
        t[int(m.group(1))] = m.group(2)
    return t
