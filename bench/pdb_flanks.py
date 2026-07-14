# How many peptide residues are actually RESOLVED beyond the 9-mer core in real pMHC-II crystals?
# Answers: "is core+-1 (11mer) too short for MHC-II structure modelling?"
# Canonical2026 (tcren-ms), class from orient_metadata.json. 2026-07-14
import gzip, json, os, collections

import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--structures", required=True, help="Canonical2026 dir of <pdb>.pdb.gz")
_ap.add_argument("--metadata", required=True, help="orient_metadata.json (has mhc.class, species)")
_a = _ap.parse_args()
D, META = _a.structures, _a.metadata
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}


def chains(path):
    """{chain: [(resseq, aa, {atom: (x,y,z)})]} for polymer residues, in file order."""
    out = collections.defaultdict(dict)
    with gzip.open(path, "rt", errors="ignore") as fh:
        for L in fh:
            if not (L.startswith("ATOM") or L.startswith("HETATM")):
                continue
            rn = L[17:20].strip()
            if rn not in AA3:
                continue
            alt = L[16]
            if alt not in (" ", "A"):
                continue
            ch, ri = L[21], L[22:27].strip()
            xyz = (float(L[30:38]), float(L[38:46]), float(L[46:54]))
            res = out[ch].setdefault(ri, [AA3[rn], {}])
            res[1][L[12:16].strip()] = xyz
    return {c: [(k, v[0], v[1]) for k, v in r.items()] for c, r in out.items()}


def d2(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b))


def analyse(pdb):
    cs = chains(os.path.join(D, pdb + ".pdb.gz"))
    if not cs:
        return None
    # peptide = shortest polymer chain of plausible ligand length
    cand = [(c, r) for c, r in cs.items() if 7 <= len(r) <= 30]
    if not cand:
        return None
    pc, pep = min(cand, key=lambda x: len(x[1]))
    # groove = all long chains (>=140 res) that contact the peptide
    groove = []
    for c, r in cs.items():
        if c == pc or len(r) < 140:
            continue
        near = sum(1 for _, _, at in r for a in at.values()
                   for _, _, pat in pep for p in pat.values() if d2(a, p) < 36)
        if near > 0:
            groove.append(r)
    if not groove:
        return None
    gat = [a for r in groove for _, _, at in r for a in at.values()]
    # burial per peptide residue: heavy-atom groove contacts within 5 A
    bur = []
    for _, _, at in pep:
        bur.append(sum(1 for a in at.values() for g in gat if d2(a, g) < 25))
    n = len(pep)
    if n < 9:
        return None
    s = max(range(n - 8), key=lambda i: sum(bur[i:i + 9]))     # the 9-mer core = most buried window
    seq = "".join(a for _, a, _ in pep)
    return dict(pdb=pdb, n=n, core_start=s, nflank=s, cflank=n - s - 9,
                seq=seq, core=seq[s:s + 9])


meta = {m["pdb.id"]: m for m in json.load(open(META))}
rows = collections.defaultdict(list)
for m in meta.values():
    if m.get("status") != "ok":
        continue
    r = analyse(m["pdb.id"])
    if r:
        rows[(m["mhc.class"], m.get("species", "?"))].append(r)

for cls in ("MHCII", "MHCI"):
    allr = [r for (c, sp), v in rows.items() if c == cls for r in v]
    if not allr:
        continue
    print(f"\n===== {cls}  (n={len(allr)} structures) =====")
    lens = collections.Counter(r["n"] for r in allr)
    print("  RESOLVED peptide length:", dict(sorted(lens.items())))
    if cls == "MHCII":
        nf = collections.Counter(r["nflank"] for r in allr)
        cf = collections.Counter(r["cflank"] for r in allr)
        print("  N-flank residues resolved beyond core:", dict(sorted(nf.items())))
        print("  C-flank residues resolved beyond core:", dict(sorted(cf.items())))
        import statistics
        print(f"  median resolved len {statistics.median([r['n'] for r in allr]):.0f}"
              f" | median N-flank {statistics.median([r['nflank'] for r in allr]):.0f}"
              f" | median C-flank {statistics.median([r['cflank'] for r in allr]):.0f}")
        ge = sum(1 for r in allr if r["nflank"] >= 2 and r["cflank"] >= 2)
        print(f"  structures resolving >=2 flanking residues on BOTH sides: "
              f"{ge}/{len(allr)} = {100*ge/len(allr):.0f}%")
        n11 = sum(1 for r in allr if r["n"] <= 11)
        print(f"  structures whose resolved peptide is <=11 residues: "
              f"{n11}/{len(allr)} = {100*n11/len(allr):.0f}%")
        for (c, sp), v in sorted(rows.items()):
            if c == cls:
                print(f"    {sp:8s} n={len(v):3d}  median len "
                      f"{statistics.median([r['n'] for r in v]):.0f}")
