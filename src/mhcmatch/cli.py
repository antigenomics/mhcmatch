"""Command-line interface for mhcmatch: ``mhcmatch <command> ...``.

Commands: ``decompose`` (no data needed), ``restriction``, ``scan``, ``logo`` (need a pmhc_data
table via ``--pmhc`` or ``$MHCMATCH_PMHC``), ``source`` (needs a proteome FASTA), and ``span``
(core -> full presented ligand; the panel is optional, and only supplies the observed-ligand tier).

**Every peptide-keyed command takes ``--peptides FILE`` and emits TSV**, because the expensive part
of almost all of them is setup that a per-peptide invocation pays again every time: the presentation
and affinity calibrators are ~5 s, the binder calibrator ~45 s, and a human-proteome length index
~70 s. Those are all cached on the :class:`~mhcmatch.Store` / :class:`~mhcmatch.Proteome` for the
life of the process, so one process over a whole list is the difference between 49 s per peptide and
a few thousand per second. A shell ``while read`` loop around the single-peptide form is the wrong
way to use this CLI and the benchmark repo measures exactly how wrong (``bench/cli/``).

``--threads`` is offered only where it does something: the C++ neighbour search (``source``,
``mimics``) releases the GIL and scales across cores. The scoring heads are small numpy products per
peptide, so threads there would buy nothing and the flag is not offered rather than being offered
and ignored.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

from . import Proteome, Store


def _add_store_opts(p):
    p.add_argument("--pmhc", help="pmhc_data TSV(.gz); else $MHCMATCH_PMHC/pmhc_<tier>.tsv.gz, "
                                  "else auto-fetched from the public HF dataset isalgo/pmhc_data")
    p.add_argument("--tier", default="full", choices=("full", "shortlist"))
    p.add_argument("--species", default="human", choices=("human", "mouse"))


def _add_batch_opts(p, what="peptide"):
    """``--peptides`` / ``--out``: the batch form of a per-``what`` command."""
    p.add_argument("--peptides", metavar="FILE",
                   help=f"run over many {what}s in ONE process: a file with one per line, or a TSV "
                        "with a `peptide` column (.gz ok, `-` = stdin). Output becomes TSV. This is "
                        "the fast path -- the per-call setup is paid once, not once per peptide")
    p.add_argument("--out", metavar="FILE", help="write TSV here instead of stdout")


def _add_thread_opt(p):
    p.add_argument("--threads", type=int, default=0, metavar="N",
                   help="worker threads for the C++ neighbour search (0 = every core)")


def _store(a):
    return Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species)


def _read_seq(arg):
    """A raw sequence, or the concatenated sequences of a FASTA file path."""
    if os.path.exists(arg):
        from .proteome import read_fasta
        seqs = read_fasta(arg)
        if seqs:
            return "".join(seqs.values())
    return arg.strip()


def _read_peptides(path, inline=()):
    """Peptides from ``path`` (one per line, or the ``peptide`` column of a TSV) plus any inline.

    ``-`` reads stdin, so this composes with a pipe. Whole-file reads on purpose: the scoring paths
    are vectorised or amortised over one setup, so handing them a whole deposit is both the fast
    path and the intended one."""
    peps = [p.strip().upper() for p in (inline or ()) if p and p.strip()]
    if not path:
        return peps
    if path == "-":
        fh, close = sys.stdin, False
    else:
        fh = (gzip.open if str(path).endswith(".gz") else open)(path, "rt")
        close = True
    try:
        first = fh.readline()
        cols = first.rstrip("\n").split("\t")
        col = cols.index("peptide") if "peptide" in cols else None
        if col is None:
            peps.append(first.strip().split("\t")[0].upper())
        for line in fh:
            line = line.rstrip("\n")
            if line:
                peps.append(line.split("\t")[col if col is not None else 0].strip().upper())
    finally:
        if close:
            fh.close()
    return [p for p in peps if p]


def _read_pairs(path, second="wt_peptide"):
    """``[(peptide, wt), ...]`` from a TSV with ``peptide`` + ``second`` columns, else columns 1-2.

    Exists so agretopicity comes out of the same pass as affinity: the mutant and its wild-type
    counterpart are one row of the caller's table, and splitting them across two runs means joining
    them back on a peptide string that is not a key."""
    op = gzip.open if str(path).endswith(".gz") else open
    fh = sys.stdin if path == "-" else op(path, "rt")
    try:
        cols = fh.readline().rstrip("\n").split("\t")
        if "peptide" in cols:
            i, j = cols.index("peptide"), (cols.index(second) if second in cols else None)
            rows = []
        else:
            i, j = 0, (1 if len(cols) > 1 else None)
            rows = [(cols[0].strip().upper(), cols[1].strip().upper() if j is not None else "")]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f and f[i].strip():
                rows.append((f[i].strip().upper(),
                             f[j].strip().upper() if j is not None and len(f) > j else ""))
        return rows
    finally:
        if path != "-":
            fh.close()


def _batch(a, positional=None):
    """The peptide list for a batch run, or ``None`` when this is a single-peptide invocation."""
    if not getattr(a, "peptides", None):
        return None
    return _read_peptides(a.peptides, [positional] if positional else ())


class _Out:
    """``--out`` if given, else stdout; closes only what it opened, and reports the row count."""

    def __init__(self, a, unit="row"):
        self.path = getattr(a, "out", None)
        self.fh = open(self.path, "w") if self.path else sys.stdout
        self.n = 0
        self.unit = unit

    def row(self, *cells):
        print("\t".join(str(c) for c in cells), file=self.fh)
        self.n += 1

    def header(self, *cells):
        print("\t".join(cells), file=self.fh)

    def close(self):
        if self.path:
            self.fh.close()
            print(f"# wrote {self.path}: {self.n:,} {self.unit}(s)", file=sys.stderr)


def cmd_decompose(a):
    store = Store()
    peps = _batch(a)
    if peps is None:
        d = store.decompose(a.peptide, cls=a.cls)
        print(f"peptide       {d.peptide}")
        print(f"anchors       {','.join(str(i + 1) for i in d.anchors)}")
        print(f"tcr_facing    {d.tcr_facing}")
        print(f"presentation  {d.presentation}")
        return
    out = _Out(a, "peptide")
    out.header("peptide", "anchors", "tcr_facing", "presentation")
    for p in peps:
        d = store.decompose(p, cls=a.cls)
        out.row(d.peptide, ",".join(str(i + 1) for i in d.anchors), d.tcr_facing, d.presentation)
    out.close()


def _resolve_panel_allele(store, name, cls):
    """Map a user-typed allele to a panel allele (exact, else prefix on punctuation-stripped name)."""
    pool = (cls,) if cls else ("mhc1", "mhc2")
    panel = {al for c in pool for al in store.alleles(c)}
    if name in panel:
        return name
    key = name.replace("*", "").replace(":", "")
    hits = sorted(a for a in panel if a.replace("*", "").replace(":", "").startswith(key))
    if hits:
        print(f"# resolved '{name}' -> '{hits[0]}'")
        return hits[0]
    print(f"# allele '{name}' not found in panel")
    return name


def cmd_restriction(a):
    store = _store(a)
    allele = _resolve_panel_allele(store, a.allele, a.cls) if a.allele else None
    peps = _batch(a)
    if peps is not None:
        out = _Out(a, "call")
        out.header("peptide", "rank", "allele", "vote", "enrichment", "score", "percent_rank",
                   "p_present", "band", "binder")
        for p in peps:
            for i, r in enumerate(store.restriction(p, cls=a.cls,
                                                    alleles=[allele] if allele else "all",
                                                    top=a.top, diffuse=a.diffuse,
                                                    calibrated=a.calibrated), 1):
                out.row(p, i, r.allele, f"{r.vote:.4g}", f"{r.enrichment:.4g}",
                        f"{r.anchor_score:.4g}" if r.anchor_score is not None else "",
                        f"{r.rank:.4g}" if a.calibrated else "",
                        f"{r.p_present:.4g}" if a.calibrated else "",
                        r.band if a.calibrated else "", "1" if r.binder else "0")
        out.close()
        return
    res = store.restriction(a.peptide, cls=a.cls, alleles=[allele] if allele else "all",
                            top=a.top, diffuse=a.diffuse, calibrated=a.calibrated)
    if not res:
        print("no presenting allele (no presentation-signature neighbours)")
        return
    diffuse = a.diffuse or a.calibrated
    hdr = f"{'allele':<18}{'vote':>7}{'enr':>7}" + ("{:>8}".format("score") if diffuse else "")
    if a.calibrated:
        hdr += f"{'%rank':>8}{'P':>7}{'band':>12}"
    print(hdr + f"{'binder':>8}")
    for r in res:
        line = f"{r.allele:<18}{r.vote:>7.2f}{r.enrichment:>7.1f}"
        if diffuse:
            line += f"{(r.anchor_score or 0.0):>8.2f}"
        if a.calibrated:
            line += f"{r.rank:>8.2f}{r.p_present:>7.2f}{r.band:>12}"
        print(line + f"{'yes' if r.binder else 'no':>8}")


def cmd_affinity(a):
    store = _store(a)
    allele = _resolve_panel_allele(store, a.allele, a.cls)
    am = store.affinity_model(a.cls)
    if a.peptides:
        # a `wt_peptide` (or second) column is read as the WT counterpart, so agretopicity comes out
        # of the same pass instead of needing a second run keyed back on the peptide
        pairs = _read_pairs(a.peptides)
        out = _Out(a, "peptide")
        out.header("peptide", "allele", "ic50_nm", "wt_peptide", "wt_ic50_nm", "amplitude", "dai")
        for p, wt in pairs:
            nm = am.predict_ic50(p, allele)
            if wt:
                out.row(p, allele, f"{nm:.6g}", wt, f"{am.predict_ic50(wt, allele):.6g}",
                        f"{am.amplitude(wt, p, allele):.6g}", f"{am.dai(wt, p, allele):.6g}")
            else:
                out.row(p, allele, f"{nm:.6g}", "", "", "", "")
        out.close()
        return
    nm = am.predict_ic50(a.peptide, allele)
    print(f"{a.peptide}  {allele}  predicted IC50 ~ {nm:,.0f} nM")
    if a.wt:
        nm_wt = am.predict_ic50(a.wt, allele)
        print(f"  WT {a.wt}: IC50 ~ {nm_wt:,.0f} nM   amplitude A=Kd_WT/Kd_MT = "
              f"{am.amplitude(a.wt, a.peptide, allele):.2f}   DAI = {am.dai(a.wt, a.peptide, allele):+.2f}")
    if a.structure:
        try:
            from .structure import StructureScorer
            sc = StructureScorer(pseudoseq=store.anchor_model(a.cls).ps)
            mj = sc.mj_energy(a.peptide, allele)
            if mj == mj:
                extra = f"   ΔΔG(WT→MT) = {sc.ddg(a.wt, a.peptide, allele):+.2f}" if a.wt else ""
                print(f"  structural MJ energy = {mj:.2f}{extra}")
            else:
                print("  (no structural template for this allele/length)")
        except ImportError as e:
            print(f"  (structure scoring unavailable: {e})")


def cmd_binder(a):
    store = _store(a)
    peps = _batch(a)
    if peps is not None:
        # The binder calibrator is the ~45 s of this command and it is cached on the store, so the
        # whole list costs one build plus a few numpy products per peptide.
        out = _Out(a, "call")
        out.header("peptide", "rank", "allele", "binder_rank", "band", "p_binder",
                   "presentation_rank", "affinity_nm", "affinity_rank")
        for p in peps:
            hits = store.binder_score(p, alleles=(a.alleles or "all"), cls=a.cls)
            for i, b in enumerate(hits[:(a.top or 10)], 1):
                out.row(p, i, b.allele, b.binder_rank, b.band, b.p_binder,
                        b.presentation_rank, b.affinity_nm, b.affinity_rank)
        out.close()
        return
    res = store.binder_score(a.peptide, alleles=(a.alleles or "all"), cls=a.cls)
    if not res:
        print("# no scorable allele (unknown groove / no background)")
        return
    print(f"# {a.peptide}: generalized binder score = geo-mean(presentation %rank, affinity %rank); "
          "lower = stronger")
    print(f"{'allele':14s}{'binder%rank':>12s}{'band':>13s}{'P(binder)':>11s}"
          f"{'pres%rank':>11s}{'aff_nM':>11s}{'aff%rank':>10s}")
    for b in res[:(a.top or 10)]:
        print(f"{b.allele:14s}{b.binder_rank:12.3f}{b.band:>13s}{b.p_binder:11.4f}"
              f"{b.presentation_rank:11.3f}{b.affinity_nm:11.0f}{b.affinity_rank:10.3f}")


def cmd_scan(a):
    hits = _store(a).scan_protein(_read_seq(a.protein), cls=a.cls or "mhc1",
                                  alleles=[a.allele] if a.allele else "all", top=a.top,
                                  correction=a.correction)
    label = f" ({a.correction} FWER/FDR)" if a.correction else ""
    print(f"# {len(hits)} presented window(s){label}")
    for pos, pep, binders in hits:
        print(f"{pos:>5}  {pep:<14}  {','.join(b.allele for b in binders)}")


def cmd_source(a):
    pm = Proteome.from_fasta(a.proteome) if os.path.exists(a.proteome) else Proteome.from_hf(a.proteome)
    peps = _batch(a)
    if peps is not None:
        # One index build per length (~70 s each for the human proteome) and one threaded C++ batch
        # query, rather than one index build per invocation.
        res = pm.find_sources(peps, max_subs=a.max_subs, exclude_exact=a.exclude_exact,
                              threads=a.threads)
        out = _Out(a, "hit")
        out.header("peptide", "protein", "position", "n_subs", "ref_peptide", "mutations")
        for p in peps:
            for h in res.get(p, ())[:(a.top or 0) or None]:
                out.row(p, h.protein, h.position, h.n_subs, h.ref_peptide,
                        ",".join(f"{q}{i + 1}{r}" for i, q, r in h.mutations) or "exact")
        out.close()
        return
    hits = pm.find_source(a.peptide, max_subs=a.max_subs, exclude_exact=a.exclude_exact)
    if not hits:
        print("# no source within max_subs")
        return
    for h in hits:
        muts = ",".join(f"{q}{i + 1}{r}" for i, q, r in h.mutations) or "exact"
        print(f"{h.protein}\tpos {h.position}\tsubs {h.n_subs}\t{h.ref_peptide}\t{muts}")


def cmd_span(a):
    from . import ligand
    prot = _read_seq(a.protein)
    corpus = None
    if a.pmhc or os.environ.get("MHCMATCH_PMHC"):
        corpus = _store(a)._panel["mhc2"].epitopes
    sp = ligand.presented_span(a.core.strip().upper(), prot, corpus=corpus, mode=a.mode,
                               flanks=tuple(int(x) for x in a.flanks.split(",")))
    if sp is None:
        print("# no reference ligand contains this core (mode=observed)")
        return
    nl, nr = sp.flanks
    print(f"tier      {sp.source}")
    print(f"core      {sp.core} @ {sp.core_start}")
    print(f"peptide   {sp.peptide}  [{sp.start}:{sp.end}]  len {len(sp.peptide)}")
    print(f"flanks    {nl} / {nr}" + (f"   clipped {sp.clipped}" if any(sp.clipped) else ""))
    print(f"score     {sp.score:+.2f}")
    print(f"alts      {sp.n_alternatives}" + (f"   support {sp.support}" if sp.support else ""))


def cmd_predict(a):
    from . import predict as P
    store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
    alleles = [x.strip() for x in a.alleles.split(",") if x.strip()]
    preds = P.predict_fasta(store, a.cls, a.fasta, alleles, rank_threshold=a.rank_threshold,
                            top=a.top, background=a.background, footprint=a.footprint, seed=a.seed)
    if a.native:
        P.write_native(preds, a.native)
        print(f"# wrote {a.native}")
    if a.scored_csv:
        P.write_scored_csv(preds, a.scored_csv)
        print(f"# wrote {a.scored_csv}")
    if not a.native and not a.scored_csv:
        print(f"# {len(preds)} predicted binder(s) (%rank <= {a.rank_threshold}) over "
              f"{len(alleles)} allele(s)")
        for p in preds[:(a.top or 20)]:
            print(f"{p.peptide:<15} {p.allele:<18} %rank={p.percent_rank:<6} {p.band:<11} "
                  f"{p.var.get('gene_name', '')}")


def cmd_logo(a):
    from . import logo
    m = logo.motif(_store(a), a.allele, a.cls or "mhc1")
    print(f"# {a.allele}  width={m['width']}  n={m['n']}  lengths={dict(sorted(m['length_hist'].items()))}")
    for i, (bits, col) in enumerate(zip(m["bits"], m["pwm"]), 1):
        top = sorted(col.items(), key=lambda x: -x[1])[:3]
        print(f"  pos {i:>2}  {bits:4.2f} bits  " + " ".join(f"{aa}:{p:.2f}" for aa, p in top))


def _read_alleles(arg):
    """``--alleles`` accepts a comma-separated list or a path to a file holding one (pipeline form)."""
    if arg and os.path.exists(arg):
        arg = open(arg).read().strip()
    return [x.strip() for x in (arg or "").replace("\n", ",").split(",") if x.strip()]


def _load_refs(spec):
    """``name=path[,name=path]`` -> ``{name: {peptides}}`` for the exact-match known-epitope flag.

    Each file is read one peptide per line, or as a TSV whose first column is the peptide."""
    refs = {}
    for part in (x.strip() for x in (spec or "").split(",") if x.strip()):
        name, _, path = part.partition("=")
        if not path:
            name, path = os.path.splitext(os.path.basename(name))[0], name
        peps = set()
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                p = line.split("\t")[0].strip().upper()
                if p and p.isalpha():
                    peps.add(p)
        refs[name] = peps
    return refs


def cmd_rank(a):
    """Rank neoantigen candidates from a window FASTA or an already-scored table."""
    from . import rank as R
    # None -> mhcmatch.known's built-in sets; --no-known-refs -> {} -> lookup off
    refs = _load_refs(getattr(a, "refs", None)) if getattr(a, "refs", None) else \
        ({} if getattr(a, "no_known_refs", False) else None)
    if a.mode == "fasta":
        store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
        rows = R.rank_fasta(store, a.input, _read_alleles(a.alleles), cls=a.cls,
                            tissue=a.tissue, tumor=a.tumor, refs=refs,
                            rank_threshold=a.rank_threshold)
    else:
        store = None
        if a.recompute_presentation:
            store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
        rows = R.rank_table(a.input, tissue=a.tissue, tumor=a.tumor, refs=refs,
                            store=store, cls=a.cls)
    rows = rows[:a.top] if a.top else rows
    cols = ["rank", "peptide", "allele", "gene", "score", "presentation", "agretopicity",
            "physchem", "expression", "expr_imputed", "wt_peptide", "known_epitope"]
    out = open(a.out, "w") if a.out else sys.stdout
    try:
        print("\t".join(cols), file=out)
        for i, r in enumerate(rows, 1):
            print("\t".join([str(i), r.peptide, r.allele, r.gene, f"{r.score:.6g}",
                             f"{r.presentation:.4g}", f"{r.agretopicity:.4g}",
                             f"{r.physchem:.4g}", f"{r.expression:.4g}",
                             "1" if r.expression_imputed else "0", r.wt_peptide,
                             r.known_epitope]), file=out)
    finally:
        if a.out:
            out.close()
            print(f"# wrote {a.out}: {len(rows)} candidate(s)")
    n_known = sum(1 for r in rows if r.known_epitope)
    if n_known:
        print(f"# {n_known} candidate(s) matched a known-epitope reference exactly "
              "(sorted into the top tier)", file=sys.stderr)


def cmd_explain(a):
    """Print every component of the aggregate for one (peptide, allele), so a rank is auditable."""
    from . import complement as CM, ipred, posbayes, rank as R
    store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
    from . import predict as P
    if a.peptides:
        pairs = _read_pairs(a.peptides)
        peps = [p for p, _ in pairs]
        recog = CM.score(peps, a.species)                 # vectorised: one pass for the whole list
        llr = [posbayes.llr(p, a.species) for p in peps]
        am = store.affinity_model(a.cls) if any(w for _, w in pairs) else None
        out = _Out(a, "peptide")
        out.header("peptide", "allele", "presentation_rank", "affinity_nm", "affinity_rank",
                   "binder_rank", "presentation_term", "recognition", "posbayes_llr", "ipred_logp",
                   "wt_peptide", "dai", "aggregate_p")
        for (p, wt), rc, lr in zip(pairs, recog, llr):
            bs = P.binder_score(store, p, alleles=[a.allele], cls=a.cls)
            b = bs[0] if bs else None
            pres = R._neglog10(b.binder_rank) if b else float("nan")
            out.row(p, a.allele,
                    f"{b.presentation_rank:.6g}" if b else "", f"{b.affinity_nm:.6g}" if b else "",
                    f"{b.affinity_rank:.6g}" if b else "", f"{b.binder_rank:.6g}" if b else "",
                    f"{pres:.6g}", f"{rc:.6g}", f"{lr:.6g}", f"{ipred.log_p(p):.6g}",
                    wt, f"{am.dai(wt, p, a.allele):.6g}" if (am and wt) else "",
                    f"{R.gate_probability(pres, rc):.6g}")
        out.close()
        return
    bs = P.binder_score(store, a.peptide, alleles=[a.allele], cls=a.cls)
    pres = R._neglog10(bs[0].binder_rank) if bs else float("nan")
    # This must be the axis `rank` scores with, and the axis GATE was fitted on -- otherwise the
    # aggregate below is a coefficient applied to a scale it was not fitted for.
    recog = R._recognition(a.peptide, a.species, cls=a.cls)
    print(f"# {a.peptide}  {a.allele}  ({a.cls})")
    if bs:
        b = bs[0]
        print(f"  presentation %rank   {b.presentation_rank:.4g}")
        print(f"  affinity   IC50 nM   {b.affinity_nm:.4g}   (%rank {b.affinity_rank:.4g})")
        print(f"  binder     %rank     {b.binder_rank:.4g}   -> presentation term {pres:+.4f}")
    print(f"  recognition  log-odds {recog:+.4f}   (complement -- the term `rank` uses)")
    print(f"  posbayes     LLR      {posbayes.llr(a.peptide, a.species):+.4f}   "
          f"(role identity only; the `aa` block of the above, shown for comparison)")
    print(f"  ipred        log P    {ipred.log_p(a.peptide):+.4f}   "
          f"(P = {ipred.p_immunogenic(a.peptide):.4f}; pooled physchem, shown for comparison)")
    if a.prior:
        # The log-odds carries no prior, so the base rate is the caller's to supply -- a screen at
        # 4.2e-4 and the training corpus at 3.2e-2 differ by ~75x.
        print(f"    P at prior {a.prior:g}   "
              f"{CM.posterior([a.peptide], a.prior, a.species)[0]:.6g}")
    if a.wt:
        from .affinity import AffinityModel
        am = AffinityModel.load(store.anchor_model(a.cls), store.corpus(a.cls), a.cls)
        print(f"  agretopicity  DAI    {am.dai(a.wt, a.peptide, a.allele):+.4f}  vs WT {a.wt}")
    if a.gene and (a.tissue or a.tumor):
        from . import expression as EX
        rec = EX.lookup(a.peptide if a.tumor else a.gene, tissue=a.tissue, tumor=a.tumor)
        print(f"  expression           {rec['median_tpm']:.4g} TPM "
              f"({'TCGA ' + a.tumor if a.tumor else 'GTEx ' + a.tissue}, n={rec['n']})"
              if rec else "  expression           (no reference row)")
    print(f"  AGGREGATE  P         {R.gate_probability(pres, recog):.6f}")


def cmd_complement(a):
    """Score peptides on the complementarity axis. Vectorised: pass the whole file, not a loop."""
    from . import complement as CM
    peps = _read_peptides(a.peptides, a.input)
    if not peps:
        raise SystemExit("no peptides: pass them as arguments or with --peptides FILE")
    s = CM.score(peps, a.species)
    out = open(a.out, "w") if a.out else sys.stdout
    try:
        if a.features:
            names = CM.feature_names(a.species)
            print("\t".join(["peptide", "score"] + names), file=out)
            X = CM.design(peps, a.species)
            for p, v, row in zip(peps, s, X):
                print("\t".join([p, f"{v:.6g}"] + [f"{x:.6g}" for x in row]), file=out)
        else:
            head = ["peptide", "score"] + (["posterior"] if a.prior else [])
            print("\t".join(head), file=out)
            post = CM.posterior(peps, a.prior, a.species) if a.prior else None
            for i, (p, v) in enumerate(zip(peps, s)):
                cells = [p, f"{v:.6g}"] + ([f"{post[i]:.6g}"] if a.prior else [])
                print("\t".join(cells), file=out)
    finally:
        if a.out:
            out.close()
            print(f"# wrote {a.out}: {len(peps)} peptide(s)", file=sys.stderr)
    if not a.prior:
        print("# score is a log-odds and carries NO prior; pass --prior to get a probability "
              f"(training corpus prevalence {CM.table(a.species)['prevalence']:.4g})",
              file=sys.stderr)


def cmd_mimics(a):
    """Near-identical reference peptides per category, in one threaded batch."""
    from . import mimics as M
    peps = _read_peptides(getattr(a, "peptides", None), a.input)
    if not peps:
        raise SystemExit("no peptides: pass them as arguments or with --peptides")
    cats = [c.strip() for c in a.categories.split(",") if c.strip()]
    unknown = [c for c in cats if c not in M.DEFAULT_REFS and c not in M.PROTEOME_REFS]
    if unknown:
        raise SystemExit(f"unknown categor(y|ies) {unknown}; expected from "
                         f"{sorted(set(M.DEFAULT_REFS) | set(M.PROTEOME_REFS))}")
    proteomes = tuple(c for c in cats if c in M.PROTEOME_REFS)
    refs = {k: v for k, v in M.DEFAULT_REFS.items() if k in cats}
    self_set, foreign = M.load_reference_sets(None, a.cls, a.species, refs=refs,
                                              proteomes=proteomes)
    ref_sets = ({"thymus": self_set} if self_set else {}) | foreign
    near = M.neighbours(peps, ref_sets, max_subs=a.max_subs, threads=a.threads)
    exact = {c: set(v) for c, v in ref_sets.items()}
    out = _Out(a, "hit")
    out.header("peptide", "category", "kind", "n_exact", "n_near", "top_mimic", "top_subs")
    for p in peps:
        got = near.get(p, {})
        for cat in ref_sets:
            hits = [h for h in got.get(cat, ()) if h[0] <= a.near_subs]
            n_ex = 1 if p in exact[cat] else 0
            if not hits and not n_ex:
                continue
            top_subs, top = (0, p) if n_ex else hits[0]
            out.row(p, cat, M.KINDS.get(cat, "?"), n_ex, len(hits), top, top_subs)
    out.close()


def cmd_expression(a):
    """Reference expression for a gene in a normal tissue, or a peptide in a tumour type."""
    from . import expression as EX
    if a.list_contexts:
        print("# GTEx tissues:")
        for t in EX.tissues():
            print(f"  {t}")
        print("# TCGA tumour types:")
        print("  " + ", ".join(EX.tumor_types()))
        return
    rec = EX.lookup(a.key, tissue=a.tissue, tumor=a.tumor)
    if not rec:
        print(f"# no reference row for {a.key!r} in "
              f"{a.tumor or a.tissue!r}")
    else:
        print(f"{a.key}\t{a.tumor or a.tissue}\t{rec['source']}\t"
              f"median {rec['median_tpm']:.4g}\tIQR {rec['q25_tpm']:.4g}-{rec['q75_tpm']:.4g}"
              f"\tn={rec['n']}")
    if a.safety:
        print(f"# {a.key} across normal tissues (highest first):")
        for t, v in EX.safety_profile(a.key, top=a.top or 10):
            print(f"  {v:10.4g}  {t}")


#: Deposits the docs, notebooks and reference lookups read, by repo-relative path. Named here so
#: ``mhcmatch bootstrap --reference`` pre-stages a container or an offline run in one call, rather
#: than each page discovering its own download on first use.
REFERENCE_FILES = (
    "immunogenicity/chowell_rebuilt.tsv.gz",      # immunogenic vs presented self
    "immunogenicity/kesmir_rebuilt.tsv.gz",       # immunogenic vs presented non-self
    "thymus/thymus_immunopeptidome.tsv.gz",       # tolerance reference for mimicry
    "ligandome/viral_foreign_iedb.tsv.gz",        # foreign reference for mimicry
    "expression/reference_expression.tsv.gz",     # GTEx tissue + TCGA tumour medians (~105 MB)
)


def cmd_bootstrap(a):
    from .store import PMHC_REPO, fetch_file, fetch_pmhc, fetch_proteome
    tiers = ("full", "shortlist") if a.tier == "all" else (a.tier,)
    for t in tiers:
        print(f"# pmhc {t}: {fetch_pmhc(t)}")
    for name in (x.strip() for x in (a.proteome or "").split(",") if x.strip()):
        print(f"# proteome {name}: {fetch_proteome(name)}")
    for rel in (REFERENCE_FILES if a.reference else ()):
        print(f"# reference {rel}: {fetch_file(rel)}")
    print(f"# cached from HF dataset {PMHC_REPO}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mhcmatch", description="peptide-MHC presentation tools")
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("decompose", help="split a peptide into anchor / TCR-facing parts (X masks)")
    d.add_argument("peptide", nargs="?", default="")
    d.add_argument("--cls", choices=("mhc1", "mhc2"))
    _add_batch_opts(d)
    d.set_defaults(fn=cmd_decompose)

    r = sub.add_parser("restriction", help="rank presenting alleles for a peptide")
    r.add_argument("peptide", nargs="?", default="")
    r.add_argument("--allele", help="restrict to a single allele")
    r.add_argument("--cls", choices=("mhc1", "mhc2"))
    r.add_argument("--diffuse", action="store_true", help="rare-allele-aware (diffusion-shrunk anchors)")
    r.add_argument("--calibrated", action="store_true",
                   help="add per-allele %%rank, P(present), and binding band (implies --diffuse)")
    r.add_argument("--top", type=int, default=10)
    _add_store_opts(r)
    _add_batch_opts(r)
    r.set_defaults(fn=cmd_restriction)

    af = sub.add_parser("affinity", help="predict IC50 (nM) + neoantigen amplitude/DAI for a peptide")
    af.add_argument("peptide", nargs="?", default="")
    af.add_argument("--allele", required=True)
    af.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    af.add_argument("--wt", help="wild-type peptide -> also report amplitude A=Kd_WT/Kd_MT and DAI")
    af.add_argument("--structure", action="store_true",
                    help="also compute the tcren MJ contact energy / ΔΔG (needs the [structure] extra)")
    _add_store_opts(af)
    _add_batch_opts(af)
    af.set_defaults(fn=cmd_affinity)

    bd = sub.add_parser("binder",
                        help="generalized binder score (presentation x affinity) ranked over alleles")
    bd.add_argument("peptide", nargs="?", default="")
    bd.add_argument("--alleles", help="comma-separated alleles (default: the whole panel)")
    bd.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    bd.add_argument("--top", type=int, default=10)
    _add_store_opts(bd)
    _add_batch_opts(bd)
    bd.set_defaults(fn=cmd_binder)

    s = sub.add_parser("scan", help="find presented peptides in a protein (sequence or FASTA path)")
    s.add_argument("protein")
    s.add_argument("--allele")
    s.add_argument("--cls", choices=("mhc1", "mhc2"))
    s.add_argument("--top", type=int, default=3)
    s.add_argument("--correction", choices=("bonferroni", "bh"),
                   help="multiple-testing control over windows x alleles (FWER / BH-FDR)")
    _add_store_opts(s)
    s.set_defaults(fn=cmd_scan)

    so = sub.add_parser("source", help="find the self peptide a neoantigen derives from")
    so.add_argument("peptide", nargs="?", default="")
    so.add_argument("--proteome", required=True,
                    help="reference proteome FASTA(.gz) path, or an HF name auto-fetched from the "
                         "public dataset (human / mouse / a pathogen stem)")
    so.add_argument("--max-subs", type=int, default=1)
    so.add_argument("--exclude-exact", action="store_true",
                    help="drop 0-mismatch hits -- the wild-type-origin question, not the "
                         "is-this-self one")
    so.add_argument("--top", type=int, help="keep only the N nearest hits per peptide")
    _add_batch_opts(so)
    _add_thread_opt(so)
    so.set_defaults(fn=cmd_source)

    lg = sub.add_parser("logo", help="motif logo (information content) + length distribution")
    lg.add_argument("allele")
    lg.add_argument("--cls", choices=("mhc1", "mhc2"))
    _add_store_opts(lg)
    lg.set_defaults(fn=cmd_logo)

    sp = sub.add_parser("span", help="extend an MHC-II binding core to the full presented ligand")
    sp.add_argument("core", help="the 9-mer binding core")
    sp.add_argument("--protein", required=True, help="source protein sequence, or a FASTA path")
    sp.add_argument("--mode", default="auto", choices=("auto", "observed", "modeled", "fixed"))
    sp.add_argument("--flanks", default="3,3", help="left,right sizes for --mode fixed")
    _add_store_opts(sp)                 # only used to supply the observed-ligand corpus
    sp.set_defaults(fn=cmd_span)

    pr = sub.add_parser("predict", help="score a variant peptide-window FASTA -> native + .scored.csv")
    pr.add_argument("fasta", help="a .peptide.fasta (pipeline schema)")
    pr.add_argument("--alleles", required=True, help="comma-separated HLA alleles (pipeline form)")
    pr.add_argument("--cls", required=True, choices=("mhc1", "mhc2"))
    pr.add_argument("--native", help="write the native TSV here")
    pr.add_argument("--scored-csv", dest="scored_csv", help="write the pipeline .scored.csv here")
    pr.add_argument("--rank-threshold", type=float, default=2.0, help="keep binders with %%rank <= this")
    pr.add_argument("--top", type=int, help="cap binders kept per window (strongest first)")
    pr.add_argument("--background", default="proteome", choices=("ligand", "proteome", "markov"))
    pr.add_argument("--footprint", default="adaptive", choices=("anchor", "core", "adaptive"))
    pr.add_argument("--seed", type=int, default=0)
    _add_store_opts(pr)
    pr.set_defaults(fn=cmd_predict)

    bs = sub.add_parser("bootstrap", help="pre-fetch the pmhc panel (and optionally proteomes) from HF")
    bs.add_argument("--tier", default="all", choices=("full", "shortlist", "all"),
                    help="which panel tier(s) to download into the huggingface_hub cache")
    bs.add_argument("--proteome", help="also fetch these reference proteomes "
                                       "(comma-separated: human,mouse,<pathogen stem>)")
    bs.add_argument("--reference", action="store_true",
                    help="also fetch the corpora and reference tables the docs, notebooks and "
                         "expression/mimicry lookups read (~115 MB) — everything offline in one call")
    bs.set_defaults(fn=cmd_bootstrap)

    rk = sub.add_parser("rank", help="rank neoantigen candidates (FASTA of windows, or a scored table)")
    rk.add_argument("mode", choices=("fasta", "table"),
                    help="fasta: mutation-spanning window FASTA + donor alleles. "
                         "table: a .scored.csv already produced by another tool")
    rk.add_argument("input", help="the .peptide.fasta or the .scored.csv")
    rk.add_argument("--alleles", help="comma-separated HLA alleles, or a file holding them "
                                      "(required for mode=fasta)")
    rk.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    rk.add_argument("--tissue", help="GTEx tissue for reference expression, e.g. 'Skin - Sun "
                                     "Exposed (Lower leg)' (the safety read)")
    rk.add_argument("--tumor", help="TCGA cancer_type for tumour expression, e.g. SKCM (melanoma)")
    rk.add_argument("--no-known-refs", action="store_true",
                    help="switch the exact-match flag off entirely (default: the built-in sets "
                         "from mhcmatch.known -- confirmed neoantigens, screened-negative "
                         "neoantigens, IEDB-immunogenic, thymic self, viral)")
    rk.add_argument("--refs", help="override the built-in known-epitope sets: "
                                   "name=path[,name=path]; one peptide per line or TSV col 1")
    rk.add_argument("--rank-threshold", type=float, default=2.0,
                    help="keep binders with presentation %%rank <= this (mode=fasta)")
    rk.add_argument("--recompute-presentation", action="store_true",
                    help="mode=table: rescore presentation with mhcmatch instead of trusting the "
                         "table's own columns")
    rk.add_argument("--top", type=int, help="print only the top N candidates")
    rk.add_argument("--out", help="write TSV here instead of stdout")
    _add_store_opts(rk)
    rk.set_defaults(fn=cmd_rank)

    ex = sub.add_parser("explain", help="every component of the aggregate for one (peptide, allele)")
    ex.add_argument("peptide")
    ex.add_argument("--allele", required=True)
    ex.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ex.add_argument("--wt", help="wild-type counterpart -> also report agretopicity (DAI)")
    ex.add_argument("--gene", help="source gene symbol, for the expression lookup")
    ex.add_argument("--tissue", help="GTEx tissue")
    ex.add_argument("--tumor", help="TCGA cancer_type, e.g. SKCM")
    ex.add_argument("--prior", type=float,
                    help="base rate for the recognition posterior (e.g. 4.8e-4 for an exome screen)")
    _add_store_opts(ex)
    _add_batch_opts(ex)
    ex.set_defaults(fn=cmd_explain)

    cm = sub.add_parser("complement",
                        help="complementarity score (recognition axis) for peptides — vectorised")
    cm.add_argument("input", nargs="*", help="peptide(s); or use --peptides")
    cm.add_argument("--peptides", help="file of peptides: one per line, or a TSV with a "
                                       "`peptide` column (.gz ok). Pass the whole deposit — "
                                       "scoring is vectorised, not a per-peptide loop")
    cm.add_argument("--prior", type=float,
                    help="base rate of the setting being scored, e.g. 4.2e-4 for the NCI screen. "
                         "Without it only the prior-free log-odds is printed, on purpose")
    cm.add_argument("--species", default="human", choices=("human", "mouse"),
                    help="which fitted table to score with; the two hosts are never pooled")
    cm.add_argument("--features", action="store_true", help="emit the full design matrix too")
    cm.add_argument("--out", help="write TSV here instead of stdout")
    cm.set_defaults(fn=cmd_complement)

    mi = sub.add_parser("mimics",
                        help="near-identical reference peptides per category (self / thymus / "
                             "viral / bacterial / neoag) -- batched and threaded")
    mi.add_argument("input", nargs="*", help="peptide(s); or use --peptides")
    mi.add_argument("--categories", default="thymus,viral,neoag",
                    help="comma-separated: thymus, viral, neoag (deposits) and self, bacterial "
                         "(reference proteomes). They answer different questions -- see "
                         "mhcmatch.mimics.KINDS -- and are never summed into one score")
    mi.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    mi.add_argument("--species", default="human", choices=("human", "mouse"))
    mi.add_argument("--max-subs", type=int, default=2, help="fuzzy search radius")
    mi.add_argument("--near-subs", type=int, default=2, help="count hits within this many subs")
    _add_batch_opts(mi)
    _add_thread_opt(mi)
    mi.set_defaults(fn=cmd_mimics)

    xp = sub.add_parser("expression", help="reference expression by normal tissue or tumour type")
    xp.add_argument("key", nargs="?", default="", help="gene symbol (with --tissue) or peptide "
                                                       "(with --tumor)")
    xp.add_argument("--tissue", help="GTEx tissue (gene-keyed, normal)")
    xp.add_argument("--tumor", help="TCGA cancer_type (peptide-keyed, tumour); SKCM = melanoma")
    xp.add_argument("--safety", action="store_true",
                    help="also print the gene across normal tissues, highest first")
    xp.add_argument("--top", type=int, help="rows for --safety (default 10)")
    xp.add_argument("--list-contexts", action="store_true",
                    help="list every GTEx tissue and TCGA tumour type available")
    xp.set_defaults(fn=cmd_expression)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
