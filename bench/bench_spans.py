#!/usr/bin/env python3
# Held-out validation of the ligand-span (flank/context) model.
#
# MHC-II task: given the binding core of a held-out eluted ligand, recover the ligand's observed
# span in its source protein. MHC-I task: no span to predict (the peptide IS the ligand), so we
# score context discrimination -- real ligands vs length-matched decoy windows from the same
# proteins.
#
#   python bench/bench_spans.py --spans /tmp/spans.tsv \
#       --proteome ~/hf/pmhc_data/proteome/human.fasta.gz ~/hf/pmhc_data/proteome/mouse.fasta.gz \
#       --cls mhc2 --seed 0
# 2026-07-14
from __future__ import annotations

import argparse
import os
import random
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mhcmatch.ligand import CORE_LEN, SpanModel                      # noqa: E402
from mhcmatch.store import _mhc2_register                            # noqa: E402
from train_spans import LEN, fit, read_proteome, _open               # noqa: E402


def read_genes(paths):
    """{accession: gene} from the FASTA ``GN=`` field. Folds isoforms and TrEMBL duplicates of the
    same gene together -- the unit a split must respect (a protein-level split still leaks across
    isoforms)."""
    genes = {}
    for path in paths:
        with _open(path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                m = re.search(r"\bGN=(\S+)", line)
                genes[acc] = m.group(1) if m else acc      # no GN -> its own group
    return genes


def load_spans(path, cls):
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            c, pep, acc, s, e, allele, pmid = line.rstrip("\n").split("\t")
            if c == cls:
                out.append((pep, acc, int(s), int(e), allele, pmid))
    return out


def split_by_gene(spans, genes, frac=0.2, seed=0):
    """Group-split on GENE, never on peptide. Nested ligand sets mean a peptide-level split puts a
    ligand in train and its +/-1 sibling in test -- that measures memorisation, not generalisation.

    A gene split alone is not enough: the same peptide sequence occurs in several genes (paralogs,
    shared domains), so it can land in both folds. Those peptides are purged from TRAIN, leaving the
    test fold intact.
    """
    rng = random.Random(seed)
    gs = sorted({genes.get(sp[1], sp[1]) for sp in spans})
    rng.shuffle(gs)
    test_g = set(gs[:max(1, int(len(gs) * frac))])
    te = [sp for sp in spans if genes.get(sp[1], sp[1]) in test_g]
    te_p = {sp[0] for sp in te}
    tr = [sp for sp in spans
          if genes.get(sp[1], sp[1]) not in test_g and sp[0] not in te_p]
    n_purged = sum(1 for sp in spans
                   if genes.get(sp[1], sp[1]) not in test_g and sp[0] in te_p)
    print(f"# purged {n_purged:,} train spans whose peptide also occurs in a test gene",
          file=sys.stderr)
    return tr, te, test_g


def truth_sets(spans):
    """{(acc, core_start): {(start, end), ...}} -- every observed span containing that core.

    Nested sets are the ground truth: a core typically has SEVERAL legitimate spans, so exact-span
    top-1 is capped well below 1.0 by ambiguity alone. We report that ceiling next to every number.
    """
    by = defaultdict(set)
    for pep, acc, s, e, _, _ in spans:
        r = _mhc2_register(pep)
        if r is None:
            continue
        by[(acc, s + r)].add((s, e))
    return by


def centered(core_start, L, plen):
    nl = (L - CORE_LEN) // 2
    s = max(0, min(core_start - nl, plen - L))
    return s, min(s + L, plen)


def evaluate(model, truth, prot, mode, modal_len):
    """mode: 'full' | 'context' (no length prior) | 'length' (modal length, centered)."""
    hit = iou_t = 0
    dn, dc, n = [], [], 0
    for (acc, cs), gt in truth.items():
        p = prot.get(acc)
        if not p:
            continue
        if mode == "length":
            ps, pe = centered(cs, modal_len, len(p))
        else:
            m = model if mode == "full" else SpanModel(
                ctx=model.ctx, lens={L: 1.0 for L in model.lens},   # flat length prior
                padbg=model.padbg, background=model.background)
            ps, pe, _, _ = m.best_span(p, cs, CORE_LEN)
        n += 1
        if (ps, pe) in gt:
            hit += 1
        # boundary error to the NEAREST member of the nested set, N and C reported separately
        best = min(gt, key=lambda g: abs(ps - g[0]) + abs(pe - g[1]))
        dn.append(abs(ps - best[0]))
        dc.append(abs(pe - best[1]))
        iou_t += max((min(pe, ge) - max(ps, gs)) / (max(pe, ge) - min(ps, gs))
                     for gs, ge in gt)
    med = lambda x: sorted(x)[len(x) // 2] if x else float("nan")    # noqa: E731
    return dict(n=n, set_recall=hit / n, iou=iou_t / n,
                dN_med=med(dn), dC_med=med(dc),
                dN_mean=sum(dn) / n, dC_mean=sum(dc) / n)


def auroc(pos, neg):
    lab = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    lab.sort(key=lambda x: x[0])
    r, i = {}, 0
    while i < len(lab):
        j = i
        while j < len(lab) and lab[j][0] == lab[i][0]:
            j += 1
        for k in range(i, j):
            r[k] = (i + j - 1) / 2 + 1
        i = j
    sp = sum(r[k] for k, (_, y) in enumerate(lab) if y == 1)
    np_, nn = len(pos), len(neg)
    return (sp - np_ * (np_ + 1) / 2) / (np_ * nn)


def bench_mhc1(model, spans, prot, rng, n=4000):
    """No span to predict -- score context discrimination: real ligands in their true source
    context vs length-matched decoy windows from the SAME proteins."""
    sample = rng.sample(spans, min(n, len(spans)))
    out = {}
    for label, flank_only in (("full 12-position", False), ("flank-only 6-position", True)):
        pos, neg = [], []
        for pep, acc, s, e, _, _ in sample:
            p = prot.get(acc)
            if not p or len(p) < len(pep) + 10:
                continue
            pos.append(model.context_score(p, s, e, flank_only=flank_only))
            d = rng.randrange(0, len(p) - len(pep))
            neg.append(model.context_score(p, d, d + len(pep), flank_only=flank_only))
        out[label] = auroc(pos, neg)
    # shuffled-context control: destroy the protein context, keep the peptide -> must fall to ~0.5
    pos, neg = [], []
    for pep, acc, s, e, _, _ in sample[:1000]:
        p = prot.get(acc)
        if not p or len(p) < len(pep) + 10:
            continue
        sh = "".join(rng.sample(p, len(p)))
        i = sh.find(pep)
        i = i if i >= 0 else rng.randrange(0, len(sh) - len(pep))
        pos.append(model.context_score(sh, i, i + len(pep), flank_only=True))
        d = rng.randrange(0, len(sh) - len(pep))
        neg.append(model.context_score(sh, d, d + len(pep), flank_only=True))
    out["shuffled-context control"] = auroc(pos, neg)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spans", required=True, help="span cache from train_spans.py --spans-out")
    ap.add_argument("--proteome", required=True, nargs="+")
    ap.add_argument("--cls", default="mhc2", choices=("mhc1", "mhc2"))
    ap.add_argument("--frac", type=float, default=0.2, help="test fraction, by gene")
    ap.add_argument("--n-test", type=int, default=4000, help="test cores to evaluate")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    prot = read_proteome(a.proteome)
    genes = read_genes(a.proteome)
    spans = load_spans(a.spans, a.cls)
    tr, te, test_g = split_by_gene(spans, genes, a.frac, a.seed)

    # --- split integrity: these must hold or every number below is meaningless ---
    tr_g = {genes.get(s[1], s[1]) for s in tr}
    assert tr_g.isdisjoint(test_g), "gene leak"
    tr_p, te_p = {s[0] for s in tr}, {s[0] for s in te}
    assert tr_p.isdisjoint(te_p), "peptide leak"
    print(f"# {a.cls}: {len(spans):,} spans | train {len(tr):,} ({len(tr_g):,} genes) "
          f"| test {len(te):,} ({len(test_g):,} genes)", file=sys.stderr)
    print(f"# split integrity OK: genes disjoint, peptides disjoint", file=sys.stderr)

    model, _ = fit(tr, prot, a.cls)                       # TRAIN FOLD ONLY
    m = SpanModel(ctx=model, lens=fit(tr, prot, a.cls)[1], padbg=0.0105)
    modal = max(m.lens, key=m.lens.get)

    if a.cls == "mhc1":
        print(f"\n## MHC-I context discrimination (AUROC vs length-matched decoys)\n")
        for k, v in bench_mhc1(m, te, prot, rng).items():
            print(f"| {k:26s} | {v:.3f} |")
        return

    truth = truth_sets(te)
    keys = list(truth)
    rng.shuffle(keys)
    truth = {k: truth[k] for k in keys[:a.n_test]}
    single = sum(1 for v in truth.values() if len(v) == 1)
    ceil = single / len(truth)
    mean_spans = sum(len(v) for v in truth.values()) / len(truth)
    print(f"# {len(truth):,} test cores | mean {mean_spans:.2f} observed spans/core "
          f"| exact-span ORACLE CEILING {ceil:.3f}", file=sys.stderr)

    rows = []
    for mode, label in (("length", f"modal length ({modal}mer), centered"),
                        ("context", "context only (no length prior)"),
                        ("full", "flank model (length + context)")):
        rows.append((label, evaluate(m, truth, prot, mode, modal)))

    # leak canary: shuffle each protein -> the context signal must vanish and the model must
    # collapse to the length/centering baseline.
    shuf = {acc: "".join(rng.sample(p, len(p))) for acc, p in
            ((k[0], prot[k[0]]) for k in truth if k[0] in prot)}
    rows.append(("leak canary (shuffled proteins)",
                 evaluate(m, {k: v for k, v in truth.items() if k[0] in shuf},
                          shuf, "full", modal)))

    print(f"\n| model | set-recall | IoU | median dN | median dC | mean dN | mean dC |")
    print(f"|---|---|---|---|---|---|---|")
    for label, r in rows:
        print(f"| {label} | {r['set_recall']:.3f} | {r['iou']:.3f} | {r['dN_med']:.0f} | "
              f"{r['dC_med']:.0f} | {r['dN_mean']:.2f} | {r['dC_mean']:.2f} |")
    print(f"\nexact-span oracle ceiling {ceil:.3f} "
          f"({mean_spans:.2f} observed spans per core -- nested sets)")


if __name__ == "__main__":
    main()
