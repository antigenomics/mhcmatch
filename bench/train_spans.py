#!/usr/bin/env python3
# Fit the ligand-span (flank/context) model from an IEDB eluted-ligand dump + reference proteomes.
# Emits src/mhcmatch/data/ligand_context.tsv, consumed by mhcmatch.ligand.load_span_model().
#
# What this models: P(observed eluted-ligand span | source protein). NOT protease cleavage -- MHC-II
# is bind-first-trim-later, so there is no strong sequence-specific endoprotease step to simulate
# (Paul 2018, PMID 30127785: a dedicated MHC-II cleavage motif reaches AUC 0.767 on ligands and has
# zero power on CD4 epitopes). The field's actual instrument is a learned flank model over eluted
# ligands -- NetMHCIIpan's -context (PMID 30446001) and MHCflurry-2.0's processing model
# (PMID 32711842). This is that.
#
#   python bench/train_spans.py --iedb ~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz \
#       --proteome ~/hf/pmhc_data/proteome/human.fasta.gz \
#       --out src/mhcmatch/data/ligand_context.tsv
# 2026-07-14
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from mhcmatch.diffusion import PROTEOME_AA_FREQ  # noqa: E402
from mhcmatch.ligand import CTX_KEYS, LIGAND_KEYS, PAD  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"

# The dump has TWO header rows and 112 columns with heavily duplicated sub-header names
# ('Name' x7, 'Starting Position' x4). csv.DictReader is last-occurrence-wins, so row["Name"]
# silently returns the ALLELE and row["Starting Position"] an empty in-vitro-process column.
# Pin every column by index and assert the (group, subheader) pair at load.
COL = {"pmid": 3, "pep": 11, "start": 15, "end": 16,
       "src_iri": 20, "parent_iri": 22, "method": 90, "mhc": 107, "cls": 111}
EXPECT = {3: ("Reference", "PMID"), 11: ("Epitope", "Name"),
          15: ("Epitope", "Starting Position"), 16: ("Epitope", "Ending Position"),
          20: ("Epitope", "Source Molecule IRI"), 22: ("Epitope", "Molecule Parent IRI"),
          90: ("Assay", "Method"), 107: ("MHC Restriction", "Name"),
          111: ("MHC Restriction", "Class")}
LEN = {"mhc2": (12, 21), "mhc1": (8, 11)}
_ACC = re.compile(r"([A-Z0-9]+)(?:\.\d+)?$")


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_proteome(paths):
    """{bare UniProt accession: sequence} from UniProt-style FASTA (``>sp|ACC|NAME``)."""
    seqs, acc, buf = {}, None, []
    for path in paths:
        with _open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if acc:
                        seqs[acc] = "".join(buf)
                    parts = line[1:].split("|")
                    acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if acc:
            seqs[acc] = "".join(buf)
            acc = None
    return seqs


def accession(row):
    """UniProt accession, version stripped. Prefer 'Molecule Parent IRI' -- 12% of rows cite NCBI
    in 'Source Molecule IRI', whose coordinates are not UniProt coordinates."""
    for k in ("parent_iri", "src_iri"):
        iri = row[COL[k]]
        if "uniprot" in iri:
            m = _ACC.search(iri.rstrip("/").split("/")[-1])
            if m:
                return m.group(1)
    return ""


def parse(iedb, cls, proteome, verbose=True):
    """Stream the dump -> [(peptide, acc, start, end, allele, pmid)], coordinates RE-DERIVED.

    IEDB's annotated Start is wrong for ~8.8% of rows and *silently* wrong for ~3.8% (the peptide
    substring-matches the protein but the coordinate points elsewhere -- signal-peptide/isoform
    numbering). So the annotated coordinate is never used as a label: the peptide is located by exact
    substring match in the same FASTA used at inference, and only unique occurrences are kept.
    """
    want = "II" if cls == "mhc2" else "I"
    lo, hi = LEN[cls]
    n = Counter()
    seen, out = set(), []
    with _open(iedb) as fh:
        r = csv.reader(fh, delimiter="\t")
        g, h = next(r), next(r)
        for i, (gg, hh) in EXPECT.items():
            assert (g[i], h[i]) == (gg, hh), f"col {i}: expected {(gg, hh)}, got {(g[i], h[i])}"
        for row in r:
            n["rows"] += 1
            if row[COL["cls"]] != want:
                continue
            n["class"] += 1
            if "mass spectrometry" not in row[COL["method"]].lower():
                continue                                   # eluted ligands only; BA boundaries are
            n["el"] += 1                                   # experimenter-chosen, i.e. noise here
            pep = row[COL["pep"]].strip().upper()
            if not (lo <= len(pep) <= hi) or set(pep) - set(AA):
                continue
            n["len_ok"] += 1
            acc = accession(row)
            if not acc:
                continue
            n["has_acc"] += 1
            if (pep, acc) in seen:
                continue                                   # dedup: heavily-published ligands
            seen.add((pep, acc))                           # must not dominate the counts
            n["unique"] += 1
            prot = proteome.get(acc)
            if prot is None:
                n["acc_absent"] += 1
                continue
            k = prot.count(pep)
            if k != 1:                                     # 0 = isoform/version drift; >1 = the
                n["occ_0" if k == 0 else "occ_multi"] += 1  # flanks are ambiguous
                continue
            s = prot.find(pep)
            out.append((pep, acc, s, s + len(pep), row[COL["mhc"]], row[COL["pmid"]]))
            n["kept"] += 1
    if verbose:
        print("# attrition", file=sys.stderr)
        for k in ("rows", "class", "el", "len_ok", "has_acc", "unique",
                  "acc_absent", "occ_0", "occ_multi", "kept"):
            print(f"#   {k:12s} {n[k]:>9,}", file=sys.stderr)
    return out, n


def context(prot, s, e):
    """The 12 terminus-relative context residues of span ``[s, e)`` in ``prot``.

    3 upstream in the protein + the ligand's own first 3 + its last 3 + 3 downstream -- the
    NetMHCIIpan ``-context`` window. Positions outside the protein are ``PAD``: a ligand ending at
    the protein's C-terminus is *evidence*, not a missing value, so it is modelled rather than
    silently dropped.
    """
    def at(i):
        return prot[i] if 0 <= i < len(prot) else PAD
    return {"flankN-3": at(s - 3), "flankN-2": at(s - 2), "flankN-1": at(s - 1),
            "ligN+1": at(s), "ligN+2": at(s + 1), "ligN+3": at(s + 2),
            "ligC-3": at(e - 3), "ligC-2": at(e - 2), "ligC-1": at(e - 1),
            "flankC+1": at(e), "flankC+2": at(e + 1), "flankC+3": at(e + 2)}


def fit(spans, proteome, cls, clamp_cys=True):
    """Laplace-smoothed position frequencies over the 12 context positions + the ligand-length prior.

    The length prior is over the **ligand length** ``P(L)``, not over core-relative flank lengths:
    defining an N-/C-flank length needs a binding core, and the allele-agnostic register is tied
    across >=2 frames on ~66% of real ligands -- so a core-relative flank-length prior would encode
    our tie-breaking rule rather than biology. N/C asymmetry is carried by the context positions,
    which are fit independently per side.
    """
    ctx = {k: Counter({a: 1 for a in AA + PAD}) for k in CTX_KEYS}   # Laplace +1
    lens = Counter({L: 1 for L in range(LEN[cls][0], LEN[cls][1] + 1)})
    for pep, acc, s, e, _, _ in spans:
        for k, r in context(proteome[acc], s, e).items():
            ctx[k][r] += 1
        lens[len(pep)] += 1

    out = {}
    for k, c in ctx.items():
        tot = sum(c.values())
        f = {a: c[a] / tot for a in c}
        if clamp_cys and k in LIGAND_KEYS:
            # Cys is depleted ~8-11x at the ligand termini but NOT in the flanks -- the flanks are
            # not in the detected peptide. That is mass-spectrometry chemistry (alkylation/missed
            # ID), not processing biology. Left in, the model would refuse every Cys-containing
            # ligand. Pin C to the proteome background (log-odds 0) and rescale the rest.
            bgC = PROTEOME_AA_FREQ.get("C", 0.023)
            rest = 1.0 - f["C"]
            if rest > 0:
                scale = (1.0 - bgC) / rest
                f = {a: (bgC if a == "C" else v * scale) for a, v in f.items()}
        out[k] = f
    tot = sum(lens.values())
    return out, {L: c / tot for L, c in lens.items()}


def allele_jsd(spans, proteome, cls, min_n=5000, top=10):
    """Sanity-check the allele-agnostic pooling: JSD(per-allele context PWM, pooled PWM).

    Processing is protease biology, not groove biology, so the flank model is pooled across alleles
    (which also unlocks the ~70% of class-II eluted-ligand rows whose restriction is a placeholder
    like 'HLA class II'). Tested rather than asserted.
    """
    import math
    by = defaultdict(list)
    for sp in spans:
        a = sp[4]
        if "*" in a and ":" in a:
            by[a].append(sp)
    big = sorted((a for a, v in by.items() if len(v) >= min_n),
                 key=lambda a: -len(by[a]))[:top]
    if not big:
        return []
    pooled, _ = fit(spans, proteome, cls)
    rows = []
    for a in big:
        per, _ = fit(by[a], proteome, cls)
        ds = []
        for k in CTX_KEYS:
            p, q = per[k], pooled[k]
            m = {r: 0.5 * (p.get(r, 0) + q.get(r, 0)) for r in set(p) | set(q)}
            kl = lambda x: sum(x.get(r, 0) * math.log2(x[r] / m[r])   # noqa: E731
                               for r in x if x.get(r, 0) > 0 and m[r] > 0)
            ds.append(0.5 * kl(p) + 0.5 * kl(q))
        rows.append((a, len(by[a]), sum(ds) / len(ds), max(ds)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iedb", required=True, help="IEDB mhc_ligand_full.tsv.gz")
    ap.add_argument("--proteome", required=True, nargs="+", help="UniProt FASTA(s), .gz ok")
    ap.add_argument("--cls", default="both", choices=("mhc1", "mhc2", "both"))
    ap.add_argument("--out", required=True, help="ligand_context.tsv")
    ap.add_argument("--spans-out", help="cache the parsed spans here (TSV) for the benchmark")
    ap.add_argument("--no-clamp-cys", action="store_true", help="keep the MS Cys artifact (debug)")
    a = ap.parse_args()

    print(f"# reading proteome(s): {', '.join(a.proteome)}", file=sys.stderr)
    prot = read_proteome(a.proteome)
    print(f"# {len(prot):,} sequences", file=sys.stderr)

    classes = ("mhc1", "mhc2") if a.cls == "both" else (a.cls,)
    rows, spans_all = [], []
    for cls in classes:
        print(f"# === {cls} ===", file=sys.stderr)
        spans, _ = parse(a.iedb, cls, prot)
        if not spans:
            print(f"# WARNING: no {cls} spans; skipping", file=sys.stderr)
            continue
        for sp in spans:
            spans_all.append((cls,) + sp)
        ctx, lens = fit(spans, prot, cls, clamp_cys=not a.no_clamp_cys)
        for k in CTX_KEYS:
            for r, v in sorted(ctx[k].items()):
                rows.append((cls, f"ctx:{k}", r, f"{v:.6f}"))
        for L, v in sorted(lens.items()):
            rows.append((cls, "len", str(L), f"{v:.6f}"))

        # known-biology check: Pro is enriched INSIDE the ligand (aminopeptidase stop signal) and
        # depleted in the flank -- the opposite of the naive "Pro at flank -2" prior.
        bgP = PROTEOME_AA_FREQ["P"]
        print(f"#   Pro  ligN+2 {ctx['ligN+2']['P'] / bgP:5.2f}x   "
              f"flankN-1 {ctx['flankN-1']['P'] / bgP:5.2f}x  (vs proteome)", file=sys.stderr)
        for row in allele_jsd(spans, prot, cls):
            print(f"#   JSD vs pooled: {row[0]:20s} n={row[1]:>7,}  mean={row[2]:.4f} "
                  f"max={row[3]:.4f}", file=sys.stderr)

    with open(a.out, "w") as fh:
        fh.write("cls\tkey\tbin\tvalue\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"# wrote {a.out} ({len(rows)} rows)", file=sys.stderr)

    if a.spans_out:
        with open(a.spans_out, "w") as fh:
            fh.write("cls\tpeptide\taccession\tstart\tend\tallele\tpmid\n")
            for sp in spans_all:
                fh.write("\t".join(str(x) for x in sp) + "\n")
        print(f"# wrote {a.spans_out} ({len(spans_all):,} spans)", file=sys.stderr)


if __name__ == "__main__":
    main()
