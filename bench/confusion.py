#!/usr/bin/env python3
"""Multi-class confusion matrix: predicted vs true MHC locus, with a non-binder class.

For each held-out positive pMHC we predict the top-scoring allele (diffused anchor log-odds) and, if
no allele clears the binder gate (score > 0), call it a *non-binder*; random corpus-AA peptides are
the true non-binders. Rows = true locus / non-binder, columns = predicted, so the diagonal is correct
locus assignment plus correct non-binder rejection, and the off-diagonal shows where presentation and
binder/non-binder calls fail. Prints the table + per-class precision/recall and writes a heatmap
(``appendix/confusion_<cls>_<species>.pdf``) and ``bench/results/confusion_<cls>_<species>.md``.

    python bench/confusion.py --pmhc-dir /path --cls mhc1 --species human
"""
import argparse
import os
import random
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import load  # noqa: E402
from tune_diffusion import _LABEL, _corpus, _locus, _split, random_peptides  # noqa: E402

from mhcmatch import Store  # noqa: E402

NB = "non-binder"
_MAIN = {"mhc1": ("HLA-A", "HLA-B", "HLA-C"), "mhc2": ("DRB", "DQB", "DPB")}


def _class_of(allele, main):
    loc = _locus(allele)
    return loc if loc in main else "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--heldout", type=float, default=0.3)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--random", type=int, default=5000, help="corpus-AA non-binders")
    ap.add_argument("--fpr", type=float, default=0.05,
                    help="calibrate the binder gate to this non-binder false-positive rate "
                         "(panel-max score quantile); 0 keeps the naive score>0 gate")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "results"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    rng = random.Random(args.seed)
    sp = {"human": "HomoSapiens", "mouse": "MusMusculus"}[args.species]
    refcount = load(os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz"), args.cls, sp)
    main_loci = _MAIN[args.cls]
    classes = list(main_loci) + ["other", NB]

    test, _ = _split(refcount, args.heldout, args.cap, rng)
    label = _LABEL[args.cls]
    train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
             for a, peps in refcount.items() for p in peps if p not in test.get(a, ())]
    store = Store.from_records(train)
    panel = store.alleles(args.cls)
    model = store.anchor_model(args.cls, h=2.0, prior_strength=10)

    def top(p):
        return max((model.score(p, a), a) for a in panel)

    # calibrate the binder gate to control the non-binder false-positive rate: cutoff = the
    # (1-fpr) quantile of the panel-max score over random peptides. A single global score>0 gate
    # does not control FPR once maximized over a large panel (multiple comparisons).
    rand_top = [top(p) for p in random_peptides(*_corpus(refcount), args.random, rng)]
    srt = sorted(s for s, _ in rand_top)
    cutoff = srt[min(len(srt) - 1, int((1 - args.fpr) * len(srt)))] if args.fpr > 0 else 0.0

    def pred(s, a):
        return NB if s <= cutoff else _class_of(a, main_loci)

    conf = {t: defaultdict(int) for t in classes}
    for a, peps in test.items():                       # held-out positives, keyed by true locus
        tl = _class_of(a, main_loci)
        for p in peps:
            s, aa = top(p)
            conf[tl][pred(s, aa)] += 1
    for s, aa in rand_top:                             # true non-binders
        conf[NB][pred(s, aa)] += 1

    # report
    gate = f"FPR-calibrated cutoff={cutoff:.2f} (target non-binder FPR {args.fpr})" \
        if args.fpr > 0 else "naive score>0 gate"
    print(f"# {args.species} {label} ({args.tier}) confusion: rows=true, cols=predicted; {gate}")
    hdr = "true\\pred".ljust(11) + "".join(c[:9].rjust(11) for c in classes) + "   recall"
    print(hdr)
    lines = [hdr]
    correct = total = 0
    for t in classes:
        row = [conf[t][p] for p in classes]
        rt = sum(row)
        rec = row[classes.index(t)] / rt if rt else float("nan")
        correct += conf[t][t]
        total += rt
        s = t.ljust(11) + "".join(str(v).rjust(11) for v in row) + f"   {rec:6.3f}"
        print(s)
        lines.append(s)
    # precision per class
    prec = {}
    for c in classes:
        col = sum(conf[t][c] for t in classes)
        prec[c] = conf[c][c] / col if col else float("nan")
    print("precision".ljust(11) + "".join(f"{prec[c]:11.3f}" for c in classes))
    acc = correct / total if total else float("nan")
    print(f"# overall accuracy = {acc:.3f}  ({correct}/{total})")

    os.makedirs(args.results, exist_ok=True)
    with open(os.path.join(args.results, f"confusion_{args.cls}_{args.species}.md"), "w") as fh:
        fh.write(f"# {args.species} {label} ({args.tier}) locus + non-binder confusion "
                 f"(diffused top-1; {gate}); accuracy {acc:.3f}\n\n```\n"
                 + "\n".join(lines) + "\n```\n")

    # row-normalized heatmap (recall on the diagonal)
    out = os.path.abspath(args.out)
    stem = f"confusion_{args.cls}_{args.species}"
    with open(os.path.join(out, stem + ".dat"), "w") as fh:
        for t in classes:                       # gnuplot 'matrix' wants row-major; row=true
            rt = sum(conf[t].values()) or 1
            fh.write(" ".join(f"{conf[t][p] / rt:.4f}" for p in classes) + "\n")
    ticks = ", ".join(f"'{c}' {i}" for i, c in enumerate(classes))
    gp = os.path.join(out, stem + ".gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 4.4in,3.6in font 'Helvetica,11'
set output '{stem}.pdf'
set title '{label} {args.species}: locus + non-binder confusion (row-normalized)'
set xlabel 'predicted'; set ylabel 'true'
set xtics ({ticks}) rotate by -30; set ytics ({ticks})
set xrange [-0.5:{len(classes) - 0.5}]; set yrange [{len(classes) - 0.5}:-0.5]
set cbrange [0:1]; set cblabel 'fraction of true class'
set palette defined (0 '#f8fafc', 1 '#1d4ed8')
unset key
plot '{stem}.dat' matrix using 1:2:3 with image
""")
    subprocess.run(["gnuplot", os.path.basename(gp)], cwd=out, check=True)
    print(f"# wrote {out}/{stem}.pdf")


if __name__ == "__main__":
    main()
