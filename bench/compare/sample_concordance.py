#!/usr/bin/env python3
"""Sample concordance: mhcmatch vs NetMHCpan/NetMHCIIpan vs the pipeline's own ``.scored.csv``.

Where ``run_compare.py`` measures *accuracy* on a synthetic binder-vs-decoy task, this measures
*agreement* on a real patient sample: the ``.peptide.fasta`` windows the Gamaleya neoantigen pipeline
already produces, scored by all three predictors over the sample's own HLA alleles. It answers "if we
swap the pipeline's MHCflurry/TLimmuno2 for mhcmatch, how differently would the calls look?".

The common axis is **presentation %rank** (mhcmatch has no nM affinity -- that is a separate model):

* mhcmatch      -- ``AnchorModel.score`` -> per-allele ``RankCalibrator.percent_rank`` (lower = stronger).
* NetMHCpan     -- ``%Rank_EL`` (lower = stronger); class I NetMHCpan-4.2b, class II NetMHCIIpan-4.3i.
* pipeline ref  -- from ``.scored.csv``: class I MHCflurry ``affinity_percentile`` (lower = stronger);
  class II TLimmuno2 ``affinity`` prediction (higher = stronger; its percentile column is empty).

Every metric is computed on a **strength** orientation (higher = stronger binder), so a positive
Spearman rho always means the tools agree, regardless of each tool's native direction.

Two views:

* **A (dense, mhcmatch vs NetMHCpan)** -- tile every window into binding-length k-mers (I: 8-11,
  II: 15), score every ``(k-mer, allele)`` with both tools. Per-allele + pooled Spearman rho of the
  two %ranks; strong-binder set overlap (Jaccard at %rank <= 2 and <= 0.5); band confusion; and
  best-allele agreement over k-mers where either tool sees a binder.
* **B (3-way, on the pipeline's own calls)** -- take the ``(epitope, best_allele)`` rows the pipeline
  emitted, score those exact pairs with mhcmatch + NetMHCpan, and report pairwise Spearman among all
  three strengths plus best-allele agreement.

Allele name-spaces are bridged by ``alleles.py`` + ``mhcmatch.pseudoseq``: class-I ``best_allele``
(``HLA-A*02:01``) -> canonical ``HLA-A02:01`` via ``normalize_allele``; class-II ``best_allele`` is
already canonical (``DRB1_1301`` / ``HLA-DPA10103-DPB10401``). mhcmatch scores panel-absent alleles
zero-shot via pseudosequence diffusion; alleles unsupported by NetMHCpan or absent from the
pseudosequence set are dropped from the head-to-head and **logged** (never silently).

Run (needs ``gawk`` on PATH for NetMHCpan's ``-xls`` writer)::

    python bench/compare/sample_concordance.py --sample TESLA1 --cls both
    python bench/compare/sample_concordance.py --sample TESLA1 --cls mhc1 --limit-windows 40  # smoke

**Privacy:** TESLA1 is public -> ``bench/results/``. Alekseech is private patient data -> its outputs
and NetMHC cache go to gitignored paths (``bench/results/private/``, ``bench/compare/_cache/``); never
commit or share Alekseech-derived peptides/rows.
"""
from __future__ import annotations

import argparse
import csv
import os
import pickle
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # sibling compare
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # bench/

import alleles as al  # noqa: E402
import metrics  # noqa: E402
import netmhc  # noqa: E402
import splits  # noqa: E402
import task  # noqa: E402

from mhcmatch import Store  # noqa: E402
from mhcmatch.calibrate import RankCalibrator, band  # noqa: E402
from mhcmatch.pseudoseq import normalize_allele, resolve_allele  # noqa: E402

_AA = set("ACDEFGHIKLMNPQRSTVWY")
_LABEL = {"mhc1": "MHCI", "mhc2": "MHCII"}
_KMER_LENS = {"mhc1": (8, 9, 10, 11), "mhc2": (15,)}   # pipeline params.mhcI_epit_len / mhcII_epit_len
_TOOL = {"mhc1": "NetMHCpan-4.2b", "mhc2": "NetMHCIIpan-4.3i"}
_STRONG, _WEAK = 0.5, 2.0   # NetMHCpan %rank binding-band thresholds (shared with mhcmatch calibrate.band)

_GAMA = "/Users/mikesh/work/academy/gamaleya/epitope_pipeline"
SAMPLES = {
    "TESLA1": {"dir": os.path.join(_GAMA, "TESLA1"), "private": False},
    "Alekseech": {"dir": os.path.join(_GAMA, "Alekseech"), "private": True},
}


# ---------------------------------------------------------------- parsing ----
def parse_peptide_fasta(path: str) -> list:
    """``[(header, sequence)]`` from a pipeline ``.peptide.fasta`` (header without the leading ``>``)."""
    out, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(buf)))
                hdr, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    return out


def tile(seq: str, lengths) -> set:
    """All standard-AA windows of each length in ``lengths`` from ``seq``."""
    seq = seq.strip().upper()
    out = set()
    for L in lengths:
        for i in range(len(seq) - L + 1):
            w = seq[i:i + L]
            if _AASET_ok(w):
                out.add(w)
    return out


def _AASET_ok(w: str) -> bool:
    return bool(w) and all(c in _AA for c in w)


# ---------------------------------------------------------------- alleles ----
def to_canonical(raw: str, cls: str) -> str:
    """Sample ``best_allele`` / typing string -> harness canonical key.

    Class I strips the ``*`` (``HLA-A*02:01`` -> ``HLA-A02:01``). Class II ``best_allele`` is already
    canonical (``DRB1_1301`` / ``HLA-DPA10103-DPB10401``) so it is passed through unchanged."""
    raw = raw.strip()
    return normalize_allele(raw) if cls == "mhc1" else raw


def _scored_file(sample: str, cls: str) -> tuple:
    """``(path, delimiter)`` of the sample's ``.epitopes.scored.{csv|tsv}`` for ``cls``.

    Format drifts by sample: TESLA1 ships comma ``.csv``, Alekseech tab ``.tsv``; same columns."""
    d = SAMPLES[sample]["dir"]
    stem = f"{sample}.{'mhcI' if cls == 'mhc1' else 'mhcII'}.epitopes.scored"
    for ext, delim in ((".csv", ","), (".tsv", "\t")):
        p = os.path.join(d, stem + ext)
        if os.path.exists(p):
            return p, delim
    raise FileNotFoundError(f"no scored file for {sample} {cls} under {d}")


def sample_alleles(sample: str, cls: str) -> list:
    """Distinct sample alleles ("guess the haplotype from output files"), canonicalized.

    TESLA1: the ``best_allele`` column of the ``.scored.{csv,tsv}`` (no dedicated typing file).
    Alekseech: the HLA-LA typing table ``*_norma.alleles.tsv``."""
    if sample == "Alekseech":
        return _alekseech_alleles(SAMPLES[sample]["dir"], cls)
    scored, delim = _scored_file(sample, cls)
    keys = {}
    with open(scored) as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            a = (row.get("best_allele") or "").strip()
            if a:
                keys[to_canonical(a, cls)] = None
    return sorted(keys)


def _alekseech_alleles(d: str, cls: str) -> list:
    """Canonical class-I/II alleles from ``Alekseech_norma.alleles.tsv`` (HLA-LA per-locus calls).

    One row per (Locus, Chromosome), single ``Allele`` column in 3-field IMGT + G-group/N form
    (``A*01:01:01G``, ``DRB4*03:01N``). Reduced to two fields; **null alleles** (``N`` suffix, not
    expressed) are dropped. DQ/DP are paired α×β across the locus (all combinations, as NetMHCIIpan
    scores them — germline phasing is not available here)."""
    from mhcmatch.pseudoseq import class2_key
    path = os.path.join(d, "Alekseech_norma.alleles.tsv")
    calls = defaultdict(set)           # locus -> {two-field allele strings, e.g. 'DRB1*03:01'}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n = _two_field(row.get("Allele") or "")
            if n:
                calls[(row.get("Locus") or "").strip()].add(n)
    keys = set()
    if cls == "mhc1":
        for loc in ("A", "B", "C"):
            for a in calls.get(loc, ()):
                keys.add(normalize_allele("HLA-" + a))
    else:
        for loc in ("DRB1", "DRB3", "DRB4", "DRB5"):        # DR: beta-only
            for b in calls.get(loc, ()):
                keys.add(class2_key("", "HLA-" + b))
        for aloc, bloc in (("DQA1", "DQB1"), ("DPA1", "DPB1")):   # DP/DQ: alpha-beta pairs
            for a in calls.get(aloc, ()):
                for b in calls.get(bloc, ()):
                    keys.add(class2_key("HLA-" + a, "HLA-" + b))
    return sorted(k for k in keys if k)


def _two_field(allele: str) -> str | None:
    """``'A*01:01:01G'`` -> ``'A*01:01'``; ``None`` for a null (``N``-suffixed, non-expressed) allele."""
    import re
    allele = allele.strip()
    if "*" not in allele:
        return None
    loc, rest = allele.split("*", 1)
    two = ":".join(rest.split(":")[:2])
    if two.endswith("N"):                      # null allele -> not presented
        return None
    two = re.sub(r"[A-Za-z]+$", "", two)       # strip G/L/Q/S expression suffix
    return f"{loc}*{two}" if two else None


def coverage(cls: str, keys, panel_alleles) -> tuple:
    """Split canonical ``keys`` into ``(both, mm_only, net_only, neither)`` by scorability.

    mhcmatch can score an allele that is **in the panel** (its own reference ligands) *or* has a
    **pseudosequence** (zero-shot via diffusion) -- the two bases are independent (some panel alleles,
    e.g. ``DRB1_1301``, lack a vendored pseudosequence yet still score from their own peptides).
    NetMHCpan only scores its shipped allele list. The head-to-head uses ``both``; the rest are
    reported, never silently dropped."""
    panel = set(panel_alleles)
    both, mm_only, net_only, neither = [], [], [], []
    for k in keys:
        mm = (k in panel) or (resolve_allele(k, cls)[0] is not None)
        net = al.emit(k, cls) is not None
        (both if (mm and net) else mm_only if mm else net_only if net else neither).append(k)
    return sorted(both), sorted(mm_only), sorted(net_only), sorted(neither)


# ---------------------------------------------------------------- scoring ----
def build_mhcmatch(pmhc_dir, cls, tier, species, background, footprint):
    """One ``Store``/``AnchorModel`` for the class + the per-allele positives map for calibration."""
    rc = splits.load_canonical(pmhc_dir, cls, species, tier)
    recs = [{"epitope": p, "mhc_a": a, "mhc_class": _LABEL[cls]}
            for a, peps in rc.items() for p in peps]
    store = Store.from_records(recs)
    model = store.anchor_model(cls, footprint=footprint, background=background)
    panel = store._panel[cls]
    pos = defaultdict(list)
    for ep, a in zip(panel.epitopes, panel.alleles):
        pos[a].append(ep)
    return model, panel, pos, task.rarity(rc)


def mhcmatch_rank(model, panel, pos, alleles, peptides, seed=0):
    """``{(allele, peptide): %rank}`` (lower = stronger). Calibrated per allele, lazily."""
    cal = RankCalibrator(model, list(alleles), panel.epitopes, n=10000, seed=seed, positives=pos)
    out = {}
    for a in alleles:
        for p in peptides:
            s = model.score(p, a)
            if s == float("-inf"):
                continue
            pr = cal.percent_rank(a, s)
            if pr == pr:                       # not nan
                out[(a, p)] = pr
    return out


def netmhc_rank(alleles, peptides, cls, cache_path=None, no_cache=False):
    """``{(allele, peptide): %Rank_EL}`` (lower = stronger) via the NetMHC wrapper, pickle-cached."""
    key = ("netmhc", cls, sorted(alleles), sorted(peptides))
    if cache_path and not no_cache and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            ck, data = pickle.load(fh)
        if ck == key:
            print(f"# reused NetMHC cache ({len(data)} scores) from {cache_path}", file=sys.stderr)
            return data
    recs = netmhc.predict({a: sorted(peptides) for a in alleles}, cls, ba=False)
    data = {k: rec["rank_el"] for k, rec in recs.items() if "rank_el" in rec}
    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "wb") as fh:
            pickle.dump((key, data), fh)
    return data


def pipeline_calls(sample, cls):
    """The pipeline's own called epitopes: ``[(epitope, canonical_allele, strength)]``.

    strength is orientation-normalized (higher = stronger): class I ``-affinity_percentile``
    (MHCflurry), class II ``+affinity`` (TLimmuno2 prediction; its percentile column is empty)."""
    scored, delim = _scored_file(sample, cls)
    out = []
    with open(scored) as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            ep = (row.get("epitope") or "").strip().upper()
            a = (row.get("best_allele") or "").strip()
            if not ep or not a or not _AASET_ok(ep):
                continue
            col = "affinity_percentile" if cls == "mhc1" else "affinity"
            v = row.get(col)
            try:
                fv = float(v)
            except (TypeError, ValueError):
                fv = None
            strength = None if fv is None else (-fv if cls == "mhc1" else fv)
            out.append((ep, to_canonical(a, cls), strength))
    return out


# --------------------------------------------------------------- analysis ----
def _band_net(rank_el):
    return band(rank_el)   # same %rank thresholds apply to NetMHCpan's %Rank_EL


def view_a(mm, nm, alleles):
    """Dense mhcmatch-vs-NetMHCpan agreement over shared ``(allele, k-mer)`` pairs."""
    shared = [k for k in mm if k in nm]
    per_allele = {}
    for a in alleles:
        pairs = [(mm[k], nm[k]) for k in shared if k[0] == a]
        if len(pairs) >= 2:
            x = [-p[0] for p in pairs]   # strength = -%rank
            y = [-p[1] for p in pairs]
            per_allele[a] = (metrics.spearman(x, y), len(pairs))
    x = [-mm[k] for k in shared]
    y = [-nm[k] for k in shared]
    pooled_rho = metrics.spearman(x, y)
    # strong-binder set overlap (peptide-allele keys)
    jac = {}
    for thr in (_STRONG, _WEAK):
        sa = {k for k in shared if mm[k] <= thr}
        sb = {k for k in shared if nm[k] <= thr}
        jac[thr] = (metrics.jaccard(sa, sb), len(sa), len(sb))
    # band confusion 3x3 (rows mhcmatch, cols NetMHC)
    labels = ("strong", "weak", "non-binder")
    conf = {r: {c: 0 for c in labels} for r in labels}
    for k in shared:
        conf[band(mm[k])][_band_net(nm[k])] += 1
    return {"n": len(shared), "pooled_rho": pooled_rho, "per_allele": per_allele,
            "jaccard": jac, "confusion": conf, "labels": labels}


def best_allele_agreement(mm, nm, peptides, alleles):
    """Fraction of k-mers where mhcmatch and NetMHCpan pick the same top allele, over k-mers where
    either tool calls a binder (min %rank <= weak threshold)."""
    agree = total = 0
    for p in peptides:
        mmv = {a: mm[(a, p)] for a in alleles if (a, p) in mm}
        nmv = {a: nm[(a, p)] for a in alleles if (a, p) in nm}
        common = set(mmv) & set(nmv)
        if len(common) < 1:
            continue
        best_mm = min(common, key=lambda a: mmv[a])
        best_nm = min(common, key=lambda a: nmv[a])
        if mmv[best_mm] <= _WEAK or nmv[best_nm] <= _WEAK:   # at least one tool sees a binder
            total += 1
            agree += (best_mm == best_nm)
    return agree, total


def view_b(sample, cls, model, panel, pos, seed, cache_path, no_cache):
    """3-way agreement on the pipeline's own called ``(epitope, best_allele)`` rows."""
    calls = pipeline_calls(sample, cls)
    alleles = sorted({a for _, a, _ in calls})
    both, *_ = coverage(cls, alleles, panel.panel)
    both = set(both)
    rows = [(ep, a, st) for ep, a, st in calls if a in both and st is not None]
    if not rows:
        return None
    peps = sorted({ep for ep, _, _ in rows})
    mm = mhcmatch_rank(model, panel, pos, both, peps, seed=seed)
    nm = netmhc_rank(sorted(both), peps, cls, cache_path=cache_path, no_cache=no_cache)
    triples = [(st, mm.get((a, ep)), nm.get((a, ep)))
               for ep, a, st in rows if (a, ep) in mm and (a, ep) in nm]
    if len(triples) < 2:
        return {"n": len(triples), "rho": {}}
    pipe = [t[0] for t in triples]                 # already strength-oriented
    mms = [-t[1] for t in triples]                 # strength = -%rank
    nms = [-t[2] for t in triples]
    return {"n": len(triples),
            "rho": {"mhcmatch~netmhc": metrics.spearman(mms, nms),
                    "mhcmatch~pipeline": metrics.spearman(mms, pipe),
                    "netmhc~pipeline": metrics.spearman(nms, pipe)}}


# ----------------------------------------------------------------- report ----
def write_report(path, sample, cls, cov, a, ba, b, params):
    both, mm_only, net_only, neither = cov
    tool = _TOOL[cls]
    L = [f"# mhcmatch vs {tool} concordance — {sample} {cls}", "",
         f"Sample **{sample}** ({'private' if SAMPLES[sample]['private'] else 'public'}), class "
         f"{cls}. Common axis: presentation %rank (lower = stronger); Spearman ρ computed on a "
         f"strength orientation so **ρ > 0 = agreement**. {params}", "",
         "## Allele coverage", "",
         f"- scored by **both** ({len(both)}): {', '.join(both) or '—'}",
         f"- mhcmatch only ({len(mm_only)}): {', '.join(mm_only) or '—'}",
         f"- {tool} only ({len(net_only)}): {', '.join(net_only) or '—'}",
         f"- neither ({len(neither)}): {', '.join(neither) or '—'}", ""]
    if a:
        L += ["## View A — dense mhcmatch vs " + tool + " (all tiled k-mers × alleles)", "",
              f"- pooled Spearman ρ = **{_f(a['pooled_rho'])}** over {a['n']:,} (k-mer, allele) pairs",
              f"- strong-binder overlap (Jaccard): %rank≤{_STRONG} = {_f(a['jaccard'][_STRONG][0])} "
              f"({a['jaccard'][_STRONG][1]} mm / {a['jaccard'][_STRONG][2]} {tool}); "
              f"%rank≤{_WEAK} = {_f(a['jaccard'][_WEAK][0])} "
              f"({a['jaccard'][_WEAK][1]} mm / {a['jaccard'][_WEAK][2]} {tool})",
              f"- best-allele agreement (k-mers where either tool binds): "
              f"**{_pct(ba[0], ba[1])}** ({ba[0]}/{ba[1]})", "",
              "### Per-allele Spearman ρ", "",
              "| allele | ρ | n k-mers |", "|---|---|---|"]
        for al_, (rho, n) in sorted(a["per_allele"].items(), key=lambda kv: -_nan(kv[1][0])):
            L.append(f"| {al_} | {_f(rho)} | {n:,} |")
        L += ["", "### Band confusion (rows mhcmatch, cols " + tool + ")", "",
              "| mhcmatch＼" + tool + " | " + " | ".join(a["labels"]) + " |",
              "|---|" + "---|" * len(a["labels"])]
        for r in a["labels"]:
            L.append(f"| {r} | " + " | ".join(str(a["confusion"][r][c]) for c in a["labels"]) + " |")
        L.append("")
    if b:
        L += ["## View B — 3-way on the pipeline's own calls", "",
              f"Over {b['n']:,} pipeline-called (epitope, best_allele) rows scored by all three:", ""]
        if b.get("rho"):
            L += ["| tool pair | Spearman ρ |", "|---|---|"]
            for k, v in b["rho"].items():
                L.append(f"| {k} | {_f(v)} |")
        L.append("")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"# wrote {path}")


def _f(x):
    return "nan" if x is None or x != x else f"{x:+.3f}" if abs(x) <= 1 else f"{x:.3f}"


def _pct(k, n):
    return "n/a" if not n else f"{100.0 * k / n:.0f}%"


def _nan(x):
    return -9 if x != x else x


# ------------------------------------------------------------------- main ----
def run_one(sample, cls, args):
    print(f"\n=== {sample} {cls} ===", file=sys.stderr)
    d = SAMPLES[sample]["dir"]
    fasta = os.path.join(d, f"{sample}.{'mhcI' if cls == 'mhc1' else 'mhcII'}.peptide.fasta")
    records = parse_peptide_fasta(fasta)
    if args.limit_windows:
        records = records[:args.limit_windows]
    kmers = set()
    for _, seq in records:
        kmers |= tile(seq, _KMER_LENS[cls])
    kmers = sorted(kmers)
    print(f"# {len(records)} windows -> {len(kmers):,} unique k-mers", file=sys.stderr)

    model, panel, pos, _rarity = build_mhcmatch(
        args.pmhc_dir, cls, args.tier, args.species, args.background, args.footprint)
    alleles = sample_alleles(sample, cls)
    cov = coverage(cls, alleles, panel.panel)
    both = cov[0]
    print(f"# alleles: both={len(both)} mm_only={len(cov[1])} net_only={len(cov[2])} "
          f"neither={len(cov[3])}", file=sys.stderr)

    a_res = ba_res = None
    if both:
        cdir = os.path.join(os.path.dirname(__file__), "_cache")
        cache = os.path.join(cdir, f"concordance_{sample}_{cls}_{args.tier}_A.pkl")
        print(f"# scoring {len(kmers):,} k-mers x {len(both)} alleles ...", file=sys.stderr)
        mm = mhcmatch_rank(model, panel, pos, both, kmers, seed=args.seed)
        nm = netmhc_rank(both, kmers, cls, cache_path=cache, no_cache=args.no_cache)
        a_res = view_a(mm, nm, both)
        ba_res = best_allele_agreement(mm, nm, kmers, both)

    bcache = os.path.join(os.path.dirname(__file__), "_cache",
                          f"concordance_{sample}_{cls}_{args.tier}_B.pkl")
    b_res = view_b(sample, cls, model, panel, pos, args.seed, bcache, args.no_cache)

    priv = SAMPLES[sample]["private"]
    outdir = os.path.join(args.out, "private") if priv else args.out
    path = os.path.join(outdir, f"concordance_{sample.lower()}_{cls}.md")
    params = (f"tier={args.tier}, background={args.background}, footprint={args.footprint}, "
              f"k-mer lengths={','.join(map(str, _KMER_LENS[cls]))}, seed={args.seed}.")
    write_report(path, sample, cls, cov, a_res, ba_res, b_res, params)
    if priv:
        print(f"# NOTE: {sample} is PRIVATE — {path} is gitignored; do not commit or share.",
              file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description="mhcmatch vs NetMHCpan sample concordance")
    ap.add_argument("--sample", default="TESLA1", help="TESLA1 | Alekseech")
    ap.add_argument("--cls", default="both", choices=("mhc1", "mhc2", "both"))
    ap.add_argument("--pmhc-dir", default=os.path.expanduser("~/hf/pmhc_data"))
    ap.add_argument("--tier", default="full", choices=("full", "shortlist"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--background", default="proteome", choices=("ligand", "proteome", "markov"),
                    help="mhcmatch log-odds null; proteome = presentation axis (NetMHCpan %%Rank_EL)")
    ap.add_argument("--footprint", default="adaptive", choices=("anchor", "core", "adaptive"))
    ap.add_argument("--limit-windows", type=int, default=0, help="cap windows for a smoke run (0=all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results"))
    args = ap.parse_args(argv)
    if args.sample not in SAMPLES:
        ap.error(f"unknown sample {args.sample!r}; known: {', '.join(SAMPLES)}")
    for cls in (("mhc1", "mhc2") if args.cls == "both" else (args.cls,)):
        run_one(args.sample, cls, args)


if __name__ == "__main__":
    main()
