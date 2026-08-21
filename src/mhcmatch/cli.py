"""Command-line interface for mhcmatch: ``mhcmatch <command> ...``.

Commands: ``decompose`` (no data needed), ``restriction``, ``scan``, ``logo`` (need a pmhc_data
table via ``--pmhc`` or ``$MHCMATCH_PMHC``), ``source`` (needs a proteome FASTA), and ``span``
(core -> full presented ligand; the panel is optional, and only supplies the observed-ligand tier).

Two commands are not peptide-keyed and take neither ``--peptides`` nor a positional peptide:
``vector`` consumes a **table of cassette units** (long windows, not minimal epitopes) and emits one
cassette, and ``deslip`` consumes a **coding sequence** in nucleotides.

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

from . import Proteome, Store, pseudoseq


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


def _add_mhc2_report(p):
    """``--mhc2-report``: how much of the class-II restriction to name. Offered on the commands where
    *we* pick the allele; a command handed one echoes back what the caller typed."""
    p.add_argument("--mhc2-report", choices=pseudoseq.REPORT_MODES, default="pair",
                   help="class-II allele granularity: pair = the full alpha-beta key, NetMHCIIpan's "
                        "own naming (default); beta = the beta chain alone; isotype = DR/DP/DQ. Use "
                        "a coarser mode to compare callers -- a class-II key leads with the beta for "
                        "DR but the alpha for DP/DQ, so matching leading genes matches two different "
                        "chains")


def _allele(a, name):
    """An allele as reported: class II reduced to ``--mhc2-report``, class I untouched."""
    return (pseudoseq.class2_report(name, getattr(a, "mhc2_report", "pair"))
            if getattr(a, "cls", None) == "mhc2" else name)


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
                out.row(p, i, _allele(a, r.allele), f"{r.vote:.4g}", f"{r.enrichment:.4g}",
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
        line = f"{_allele(a, r.allele):<18}{r.vote:>7.2f}{r.enrichment:>7.1f}"
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
                out.row(p, i, _allele(a, b.allele), b.binder_rank, b.band, b.p_binder,
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
        print(f"{_allele(a, b.allele):14s}{b.binder_rank:12.3f}{b.band:>13s}{b.p_binder:11.4f}"
              f"{b.presentation_rank:11.3f}{b.affinity_nm:11.0f}{b.affinity_rank:10.3f}")


def cmd_scan(a):
    hits = _store(a).scan_protein(_read_seq(a.protein), cls=a.cls or "mhc1",
                                  alleles=[a.allele] if a.allele else "all", top=a.top,
                                  correction=a.correction)
    label = f" ({a.correction} FWER/FDR)" if a.correction else ""
    print(f"# {len(hits)} presented window(s){label}")
    for pos, pep, binders in hits:
        print(f"{pos:>5}  {pep:<14}  {','.join(_allele(a, b.allele) for b in binders)}")


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
        P.write_native(preds, a.native, core=a.core)
        print(f"# wrote {a.native}")
    if a.scored_csv:
        P.write_scored_csv(preds, a.scored_csv, core=a.core)
        print(f"# wrote {a.scored_csv}")
    if not a.native and not a.scored_csv:
        print(f"# {len(preds)} predicted binder(s) (%rank <= {a.rank_threshold}) over "
              f"{len(alleles)} allele(s)")
        for p in preds[:(a.top or 20)]:
            print(f"{p.peptide:<15} {_allele(a, p.allele):<18} %rank={p.percent_rank:<6} {p.band:<11} "
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


#: The six (component, channel) pairs of the mimicry aggregate, in artifact order.
from .rank import MIMICRY_PAIRS as _MIM_PAIRS      # noqa: E402  (the emitted-column order)


def _mimicry_scores(peptides, cls: str, no_self: bool):
    """Score a candidate list on the mimicry aggregate, warning about what it is about to cost.

    Shared by ``rank --extended/--annotate`` and ``mimicry``: the reference index is built once for
    the whole list, which is the only way this is affordable at all.

    **Only the lengths this list actually has.** The index is per-length and so is the cost -- one
    proteome pass is ~11 s, plus ~1.0 min per class-II length to resolve a register for each of
    12,685,964 windows -- and the class-II ladder is fifteen lengths deep. Indexing all of them for
    a list of 15-mers is ~19 min of work for one length's worth of answers."""
    from . import mimicry as MM
    lens = sorted({len(p) for p in peptides})
    if not no_self:
        print(f"# indexing the host proteome (~12 M windows) at {len(lens)} length(s) "
              f"{','.join(str(L) for L in lens)}: minutes and ~7 GB, paid once for the whole list. "
              "--no-self skips it at the cost of the largest coefficients",
              file=sys.stderr, flush=True)
    refs = MM.load_references(cls=cls, with_self=not no_self, lengths=lens)
    return MM.score(peptides, refs, cls=cls, allow_missing=no_self)


def _aggregate_channels(cls: str, no_self: bool, species: str = "human"):
    """Build the ``channels`` callable ``rank`` needs to score with the fitted aggregate.

    Returns ``list[peptide] -> {C_corpus_thymus, C_corpus_self, C_corpus_viral}``.

    **Since 0.24.0 there is no search here at all.** ``BOECRT`` needed ``self_tcr`` -- its
    second-largest coefficient -- so an aggregate score forced the host-proteome reference index:
    6 min 15 s and ~7.5 GB, the largest single cost in the package. 0.24.0 replaced the neighbour
    search with a :func:`mhcmatch.mimicry.corpus_spectrum` table contraction, so all three channels
    together cost three tables of 64 KB and the ranking path builds no trie at all. That is why
    ``self`` and ``viral`` are back in the model: they were dropped in 0.21.0 for what they cost,
    not for what they were worth. ``--no-self`` still matters for ``--extended``/``--annotate``,
    which report *which* reference peptide was nearest and do need the index, and for the safety
    scan.

    ``species`` picks the ``self`` proteome, so a mouse run scores mouse self. The ``thymus`` and
    ``viral`` deposits are human-only; that is stated in :func:`mhcmatch.mimicry.corpus_spectrum`
    and is a roadmap item, not a silent substitution.

    The ``C_phys`` pair is deliberately absent: :func:`mhcmatch.rank._finish` computes both, because
    they are matrix products against published residue vectors and need no deposit at all.
    """
    from . import mimicry as MM
    from . import rank as R

    def channels(peptides):
        spec = MM.corpus_spectrum(cls=cls, components=("thymus", "self", "viral"),
                                  self_species=species)
        rows = MM.corpus_R(list(peptides), spec, cls=cls)
        return {f"C_corpus_{c}": [r.get(c, float("nan")) for r in rows]
                for c in ("thymus", "self", "viral")
                if f"C_corpus_{c}" in R.CHANNEL_COLUMNS}

    return channels


def cmd_rank(a):
    """Rank neoantigen candidates from a window FASTA or an already-scored table.

    With ``--score aggregate`` (the default) every one of the model's nine features is computed
    *before* scoring and emitted as a column -- see :func:`_aggregate_channels`. Before 0.20.0 four
    of the then nine were computed after scoring, or not at all, and contributed a constant to
    every candidate while the output still said ``BOECRT``.

    ``--no-self`` and ``--score aggregate`` were mutually exclusive until 0.21.0, because
    ``BOECRT`` scored on ``self_tcr``. ``GRAND`` does not, so the combination is now allowed and
    the host-proteome index is off the ranking path entirely.
    """
    from . import rank as R
    # None -> mhcmatch.known's built-in sets; --no-known-refs -> {} -> lookup off
    refs = _load_refs(getattr(a, "refs", None)) if getattr(a, "refs", None) else \
        ({} if getattr(a, "no_known_refs", False) else None)
    if a.mode == "fasta":
        store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
        rows = R.rank_fasta(store, a.input, _read_alleles(a.alleles), cls=a.cls,
                            tissue=a.tissue, tumor=a.tumor, refs=refs,
                            rank_threshold=a.rank_threshold, score=a.score,
                            prevalence=a.prevalence,
                            channels=_aggregate_channels(a.cls, a.no_self, a.species)
                            if a.score == "aggregate" else None)
    else:
        store = None
        if a.recompute_presentation:
            store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
        rows = R.rank_table(a.input, tissue=a.tissue, tumor=a.tumor, refs=refs,
                            store=store, cls=a.cls, score=a.score,
                            prevalence=a.prevalence,
                            channels=_aggregate_channels(a.cls, a.no_self, a.species)
                            if a.score == "aggregate" else None)
    rows = rows[:a.top] if a.top else rows
    cols = list(R.BASE_COLUMNS)
    if a.score == "aggregate":
        cols += list(R.AGGREGATE_COLUMNS)
        model = rows[0].components.get("model", "") if rows else ""
        print(f"# scored with {model or 'aggregate'}: "
              f"{', '.join(R.AGGREGATE_FEATURES)}", file=sys.stderr)
    # The mimicry columns are appended, never folded into `score`. Whether mimicry belongs inside
    # the gate is a benchmark question that is not settled, and quietly moving the ranking on an
    # unvalidated term is the failure mode worth avoiding -- so the ordering is identical with and
    # without these flags.
    mim, ann = [], []
    if a.extended or a.annotate:
        mim = _mimicry_scores([r.peptide for r in rows], a.cls, a.no_self)
    if a.extended:
        cols += list(R.EXTENDED_COLUMNS)
    if a.annotate:
        from . import mimicry as MM
        ann = MM.annotate([r.peptide for r in rows], cls=a.cls)
        cols += list(R.ANNOTATE_COLUMNS)
    if a.core:
        cols += list(R.CORE_COLUMNS)
    out = open(a.out, "w") if a.out else sys.stdout
    try:
        print("\t".join(cols), file=out)
        for i, r in enumerate(rows, 1):
            cells = [str(r.rank), r.peptide, _allele(a, r.allele), r.gene, f"{r.score:.6g}",
                     f"{r.p_response:.4g}",
                     f"{r.presentation:.4g}", f"{r.binder:.4g}", f"{r.occupancy:.4g}",
                     f"{r.agretopicity:.4g}",
                     f"{r.physchem:.4g}", f"{r.expression:.4g}",
                     "1" if r.expression_imputed else "0",
                     str(r.n_alleles_presenting), r.alleles_presenting,
                     r.imputed, r.wt_peptide,
                     r.known_epitope, r.variant_type]
            if a.score == "aggregate":
                cells += [f"{r.components[c]:.6g}" for c in R.AGGREGATE_COLUMNS]
            if a.extended:
                s = mim[i - 1]
                cells += [f"{s.logodds:.6g}", f"{s.autoimmune:.6g}"]
                cells += [f"{s.components[c][ch]:.6g}" for c, ch in _MIM_PAIRS]
            if a.annotate:
                near = mim[i - 1].nearest
                for c, ch in _MIM_PAIRS:
                    n = (near.get(c) or {}).get(ch)
                    cells += [n["peptide"], n["source"], str(n["subs"])] if n else ["", "", ""]
                g = ann[i - 1]
                cells += [str(g["neoag_distance"]), g["neoag_nearest"] or "",
                          str(g["neoag_n_within"])]
            if a.core:
                cells += [r.core, str(r.core_offset), r.core_source]
            print("\t".join(cells), file=out)
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
    from . import complement as CM, posbayes, rank as R
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
                   "binder_rank", "presentation_term", "recognition", "posbayes_llr",
                   "wt_peptide", "dai", "aggregate_p")
        for (p, wt), rc, lr in zip(pairs, recog, llr):
            bs = P.binder_score(store, p, alleles=[a.allele], cls=a.cls)
            b = bs[0] if bs else None
            pres = R._neglog10(b.binder_rank) if b else float("nan")
            out.row(p, a.allele,
                    f"{b.presentation_rank:.6g}" if b else "", f"{b.affinity_nm:.6g}" if b else "",
                    f"{b.affinity_rank:.6g}" if b else "", f"{b.binder_rank:.6g}" if b else "",
                    f"{pres:.6g}", f"{rc:.6g}", f"{lr:.6g}",
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


def _read_table(path):
    """Every row of a TSV with a ``peptide`` column, as dicts, preserving column order."""
    fh = sys.stdin if path == "-" else (
        gzip.open if str(path).endswith(".gz") else open)(path, "rt")
    try:
        cols = fh.readline().rstrip("\n").split("\t")
        if "peptide" not in cols:
            raise SystemExit(f"{path}: no `peptide` column (found {cols})")
        out = []
        for line in fh:
            line = line.rstrip("\n")
            if line:
                d = dict(zip(cols, line.split("\t")))
                d["peptide"] = (d.get("peptide") or "").strip().upper()
                out.append(d)
        return out
    finally:
        if path != "-":
            fh.close()


def cmd_neoag(a):
    """Annotate candidates by fuzzy search against the tested-neoantigen database.

    **Prior evidence, not a prediction.** The database is the union of what has been assayed
    somewhere, so a hit says "this, or something one or two substitutions from it, was tested and
    came back immunogenic". That is a strong prioritisation signal for a fresh cohort and a
    meaningless one for a benchmark assembled from the same deposits, which is why
    :func:`mhcmatch.mimicry.annotate` is kept out of the fitted mimicry aggregate entirely.

    With ``--peptides`` pointing at a TSV, every original column is carried through and the
    annotation is appended, so this drops into an existing candidate table without a join."""
    from . import mimicry as MM
    src = getattr(a, "peptides", None)
    rows = _read_table(src) if src else None
    peps = [r["peptide"] for r in rows] if rows else _read_peptides(None, a.input)
    if not peps:
        raise SystemExit("no peptides: pass them as arguments or with --peptides")
    ann = MM.annotate(peps, cls=a.cls, max_subs=a.max_subs)
    cols = MM.NEOAG_COLUMNS
    extra = [c for c in (rows[0] if rows else {}) if c != "peptide"]
    # This command has no store and no allele, so a class-II core here is the allele-agnostic
    # register. `core_source` says `heuristic` rather than the docs saying it somewhere else.
    from .rank import CORE_COLUMNS
    from .store import binding_core
    out = _Out(a, "candidate")
    out.header("peptide", *extra, *cols, *(CORE_COLUMNS if a.core else ()))
    for i, r in enumerate(ann):
        if a.known_only and not r["known"]:
            continue
        if a.hits_only and r["neoag_distance"] > a.max_subs:
            continue
        s = rows[i] if rows else {}
        core = binding_core(r["peptide"], a.cls) if a.core else None
        out.row(r["peptide"], *(s.get(c, "") for c in extra),
                *("" if r[k] is None else r[k] for k in cols),
                *((core[0], core[1],
                   ("footprint" if a.cls != "mhc2" else "heuristic") if core[0] else "")
                  if a.core else ()))
    out.close()


def cmd_mimicry(a):
    """The fitted mimicry aggregate: six signed contributions, their sum, and what was hit.

    One column per (component, channel), so the aggregate can always be taken apart -- a candidate
    at +1.5 because it looks viral and one at +1.5 because it looks nothing like self are different
    candidates. ``--coefficients`` prints the shipped model instead of scoring anything.

    ``self`` is the expensive reference (~12 M windows per length); ``--no-self`` drops it, and with
    it the largest coefficients, so the aggregate there is deliberately a smaller model."""
    from . import mimicry as MM
    if a.coefficients:
        p = MM.params(a.cls)
        print(f"# mhcmatch.mimicry {a.cls}, model version {p['version']}, radius {p['radius']}")
        print(f"# fitted on {p['fit']['n']:,} rows / {p['fit']['pos']:,} positives over "
              f"{len(p['fit']['screens'])} screens, screen indicators as nuisance columns")
        print("\t".join(("component", "channel", "coef", "sd", "z")))
        for f, c, s in zip(p["features"], p["logistic"]["coef"], p["logistic"]["sd"]):
            comp, ch = f.rsplit("_", 1)
            print(f"{comp}\t{ch}\t{c:+.4f}\t{s:.4f}\t{c / s:+.2f}")
        print(f"# AUROC {p['fit']['auroc_pooled']:.3f} pooled, "
              f"{p['fit']['auroc_within_screen_median']:.3f} median within screen -- report the "
              f"second one", file=sys.stderr)
        return
    src = getattr(a, "peptides", None)
    rows = _read_table(src) if src else None
    peps = [r["peptide"] for r in rows] if rows else _read_peptides(None, a.input)
    if not peps:
        raise SystemExit("no peptides: pass them as arguments or with --peptides")
    scores = _mimicry_scores(peps, a.cls, a.no_self)
    prob = MM.probability(scores, corpus=a.corpus, cls=a.cls) if a.corpus else None
    cols = [f"{c}_{ch}" for c, ch in _MIM_PAIRS]
    if a.annotate:
        cols += [f"{k}_{c}_{ch}" for c, ch in _MIM_PAIRS for k in ("nearest", "source", "subs")]
    extra = [c for c in (rows[0] if rows else {}) if c != "peptide"]
    out = _Out(a, "peptide")
    out.header("peptide", *extra, "logodds", "autoimmune",
               *(["p_response"] if prob else []), *cols)
    for i, s in enumerate(scores):
        d = s.as_dict()
        out.row(s.peptide, *((rows[i] if rows else {}).get(c, "") for c in extra),
                f"{s.logodds:.6g}", f"{s.autoimmune:.6g}",
                *([f"{prob[i]:.6g}"] if prob else []),
                *(f"{d[c]:.6g}" if isinstance(d.get(c), float) else d.get(c, "") for c in cols))
    out.close()
    if not prob:
        print("# logodds is prior-free; --corpus screens maps it to a probability against that "
              "corpus's own prevalence, which is not transferable", file=sys.stderr)


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
        print("# TCGA tumour type -> its matched normal GTEx tissue(s), best match first.")
        print("# Keys are TCGA study abbreviations (NCI GDC); values are GTEx SMTSD names.")
        print("# Neither is a clinical coding system -- not ICD-O-3, SNOMED CT or OncoTree.")
        for t in EX.tumor_types():
            m = EX.matched_tissues(t)
            flag = "  (approximate: GTEx has no matching organ)" \
                if t in EX.TUMOR_TISSUE_APPROXIMATE else ""
            print(f"  {t:6s} {' | '.join(m) if m else '(no matched normal)'}{flag}")
        print("# CRC is not a TCGA code: TCGA has COAD and READ, and the source table merged them.")
        unmatched = [t for t in EX.tissues()
                     if not any(t in v for v in EX.TUMOR_TISSUE.values())]
        print(f"\n# {len(EX.tissues())} GTEx tissues in total; the {len(unmatched)} below are not "
              "any tumour type's matched normal and are for the safety read only:")
        for t in unmatched:
            print(f"  {t}")
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
    "immunogenicity/chowell_iedb_full.tsv.gz",    # the rebuilt corpus recognition is fitted on
    "immunogenicity/chowell_rebuilt.tsv.gz",      # immunogenic vs presented self (legacy)
    "immunogenicity/kesmir_rebuilt.tsv.gz",       # immunogenic vs presented non-self (legacy)
    "thymus/thymus_immunopeptidome.tsv.gz",       # tolerance reference for mimicry
    "ligandome/viral_foreign_iedb.tsv.gz",        # foreign reference for mimicry
    "expression/reference_expression.tsv.gz",     # GTEx tissue + TCGA tumour medians (~105 MB)
)


def _parse_quota(spec: str) -> dict:
    """``'mhc1=8:2,mhc2=4:1'`` -> ``{"mhc1": (8, 2), ...}``. Slots first, response target second."""
    out = {}
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            arm, rest = part.split("=", 1)
            slots, target = rest.split(":", 1)
            out[arm.strip()] = (int(slots), int(target))
        except ValueError:
            raise SystemExit(f"--quota: cannot read {part!r}; use ARM=SLOTS:TARGET, e.g. "
                             "'mhc1=8:2,mhc2=4:1,nonconventional=3:1'") from None
        if out[arm.strip()][1] > out[arm.strip()][0]:
            raise SystemExit(f"--quota {part}: asking for {target} responses out of {slots} slots")
    if not out:
        raise SystemExit("--quota was empty; use ARM=SLOTS:TARGET, e.g. 'mhc1=8:2,mhc2=4:1'")
    return out


def _read_units(path):
    """``[Unit]`` from a TSV with ``peptide``/``gene``/``allele``/``p`` (+ optional
    ``mutation_index``, ``cls``, ``kind``).

    ``kind`` is the variant class -- ``missense`` (the default) or a non-conventional product
    (``frameshift``, ``fusion``, ``splice``, ``retained_intron``, ``ORF``, ``editing``). ``rank``
    emits it as ``variant_type`` where the input carried one; it only matters under ``--quota``.

    Deliberately not the `rank` table read directly. `rank` emits **minimal epitopes**, and a unit is
    the long peptide around one mutation -- injecting a 9-mer would build the tolerising
    configuration (see :func:`mhcmatch.vector.unit`). The caller joins `rank`'s output back to the
    variant's protein context, which is the step that knows where the mutation sits, and this reads
    the result.
    """
    from .vector import Unit

    op = gzip.open if str(path).endswith(".gz") else open
    fh = sys.stdin if path == "-" else op(path, "rt")
    try:
        cols = fh.readline().rstrip("\n").split("\t")
        need = ("peptide", "gene", "allele", "p")
        missing = [c for c in need if c not in cols]
        if missing:
            raise SystemExit(f"{path}: missing column(s) {', '.join(missing)}; a unit table needs "
                             f"{', '.join(need)} (+ optional mutation_index, cls). `rank` gives you "
                             "gene, allele and a score -- peptide must be the long window around the "
                             "mutation, not the minimal epitope")
        ix = {c: cols.index(c) for c in cols}
        units = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[ix["peptide"]].strip():
                continue
            pep = f[ix["peptide"]].strip().upper()
            mi = (int(f[ix["mutation_index"]]) if "mutation_index" in ix
                  and len(f) > ix["mutation_index"] and f[ix["mutation_index"]].strip()
                  else len(pep) // 2)
            units.append(Unit(peptide=pep, mutation_index=mi, gene=f[ix["gene"]].strip(),
                              allele=f[ix["allele"]].strip(), p=float(f[ix["p"]]),
                              cls=(f[ix["cls"]].strip() if "cls" in ix and len(f) > ix["cls"]
                                   else "mhc1"),
                              kind=(f[ix["kind"]].strip() or "missense"
                                    if "kind" in ix and len(f) > ix["kind"] else "missense")))
        return units
    finally:
        if path != "-":
            fh.close()


def _read_unit_rows(path):
    """``[dict]`` from the same TSV :func:`_read_units` reads, without the long-window contract.

    Used only with ``--context``, where ``peptide`` is deliberately a *minimal* epitope and the long
    window is rebuilt from the FASTA instead -- so the check that belongs here is that the columns
    exist, not that the peptide is long.
    """
    op = gzip.open if str(path).endswith(".gz") else open
    fh = sys.stdin if path == "-" else op(path, "rt")
    try:
        cols = fh.readline().rstrip("\n").split("\t")
        need = ("peptide", "gene", "allele", "p")
        missing = [c for c in need if c not in cols]
        if missing:
            raise SystemExit(f"{path}: missing column(s) {', '.join(missing)}; with --context a "
                             f"candidate table needs {', '.join(need)}, which is what `rank` emits "
                             "(rename its `score` column to `p`)")
        rows = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(cols) or not f[cols.index("peptide")].strip():
                continue
            rows.append(dict(zip(cols, f)))
        return rows
    finally:
        if path != "-":
            fh.close()


def cmd_vector(a):
    """Screen, select, order: a polyepitope cassette from a table of candidate units."""
    from . import vector as V

    # Checked before anything expensive: `order` raises on this too, but only after the panel and
    # the calibrators have been built, which is ~10 s spent to learn about a typo.
    if a.objective == "rate" and a.binder_threshold is None:
        raise SystemExit("--objective rate needs --binder-threshold: the rate is binders per "
                         "register, so something has to say what counts as a binder")
    # One register vocabulary for the whole command. A class-II core is 9 residues read out of a
    # longer span, so screening it at class-I lengths would look at windows no MHC-II ever presents.
    lengths = V.JUNCTION_LENGTHS if a.cls == "mhc1" else V.MHC2_JUNCTION_LENGTHS
    if a.context:
        from .predict import parse_fasta
        rows = _read_unit_rows(a.candidates)
        records = parse_fasta(a.context)
        units = V.units_from_context(rows, records, length=a.unit_length, cls=a.cls)
        print(f"# --context: {len(rows)} ranked row(s) over {len(records)} window(s) -> "
              f"{len(units)} unit(s), one per variant", file=sys.stderr)
    else:
        units = _read_units(a.candidates)
    print(f"# {len(units)} candidate unit(s) over "
          f"{len({u.allele for u in units})} allotype(s)", file=sys.stderr)

    rejected = []
    if a.screen:
        from .proteome import gene_symbols
        from .store import fetch_proteome
        print(f"# screening: one whole-proteome index per register length ({len(lengths)} for "
              f"{a.cls}), ~12 GB peak each and a few minutes apiece. Paid once for the whole "
              "candidate list, so screen everything in one call", file=sys.stderr, flush=True)
        fa = fetch_proteome(a.species)
        risk = V.self_origin_risk(Proteome.from_fasta(fa), gene_symbols(fa, key="accession"),
                                  min_tpm=a.min_tpm, max_subs=a.max_subs)
        units, rejected = V.screen(units, risk, lengths=lengths)
        print(f"# withdrawn: {len({id(u) for u, _, _ in rejected})} unit(s), "
              f"{len(rejected)} reason(s); {len(units)} remain", file=sys.stderr)

    sel = V.select(units, n0=a.n0, cls=a.cls if a.cls_filter else None)
    print(f"# selected {len(sel.units)} of {len(units)}, expected yield "
          f"{sel.expected_yield:.2f} at n0={a.n0}", file=sys.stderr)

    comp = topk = None
    if a.quota:
        from . import portfolio as PF
        quotas = _parse_quota(a.quota)
        universe = _read_alleles(a.alleles) or sorted({u.allele for u in units if u.allele})
        comp = PF.compose(units, quotas, a.block_live, weight_evenness=a.evenness,
                          universe=universe)
        # The same slot budgets filled by score alone -- the cassette a ranked list gives you.
        # Reported beside the composed one because "different from ranking" is a claim that has to
        # be shown on the caller's own candidates, not asserted from a docstring.
        topk = {}
        for arm, (slots, target) in quotas.items():
            pool = sorted([u for u in units if PF.default_arm(u) == arm], key=lambda u: -u.p)
            topk[arm] = pool[:slots]
        print(f"# composed {len(comp.units)} unit(s) over {len(quotas)} arm(s); joint "
              f"P(all quotas met) {comp.joint:.4f} at q={a.block_live}", file=sys.stderr)

    # Which unit set becomes a *sequence*. Without `--quota` it is `select`'s, as it always was.
    # With `--quota` the composed set is the deliverable and the same slot budgets filled by score
    # alone ride along as a second cassette -- because "a portfolio is not a ranking" is a claim that
    # has to be laid out on the caller's own candidates, not asserted. Until 0.24.1 `--quota`
    # composed a set and then built the sequence from `select` anyway, so it reported and did not act.
    plans = [("cassette", sel.units)]
    if comp is not None:
        plans = [("cassette_composed", comp.units),
                 ("cassette_topk", [u for arm in quotas for u in topk[arm]])]

    built, binder, alleles = [], None, None
    for name, us in plans:
        if len(us) >= 2:
            if binder is None:
                store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
                alleles = _read_alleles(a.alleles) or sorted(
                    {u.allele for _, p_ in plans for u in p_ if u.allele})
                binder = V.store_binder(store, alleles, cls=a.cls)   # ~10 s, paid once for both
            built.append((name, V.order(us, binder=binder, lengths=lengths, alleles=alleles,
                                        objective=a.objective,
                                        binder_threshold=a.binder_threshold,
                                        threshold=a.threshold)))
        elif us:
            # One unit has no junctions, so `order` returns before it ever calls `binder` -- building
            # the panel and its calibrators here would be ~10 s spent to lay out a cassette of one.
            built.append((name, V.order(us, binder=None)))
    cas = built[0][1] if built else None

    o = _Out(a, "row")
    try:
        o.header("section", "i", "key", "value", "detail")
        # One row per (unit, register, source gene), carrying that gene's worst tissue. `screen`
        # returns a reason per tissue, and a gene is transcribed in many -- on one 7-unit test that
        # was 2,121 rows for 4 withdrawals, which is a table nobody reads. The full set is still
        # what the library returns; this is the presentation.
        worst = {}
        for u, reg, why in rejected:
            k = (id(u), reg, why.get("gene", ""))
            if k not in worst or why.get("tpm", 0) > worst[k][2].get("tpm", 0):
                worst[k] = (u, reg, why)
        for u, reg, why in worst.values():
            sub = f" {why['subs']}sub" if "subs" in why else ""
            o.row("withdrawn", "", u.gene, why.get("clause", ""),
                  f"{why.get('gene', '')}{sub} {why.get('tissue', '')} "
                  f"{why.get('tpm', 0):.1f}".strip() + (f" via {reg}" if reg else ""))
        if comp is not None:
            from . import portfolio as PF
            for arm, d in comp.arms.items():
                tk = topk[arm]
                p_top = (PF.p_at_least([u.p for u in tk], [u.allele for u in tk], a.block_live,
                                       k=d["target"]) if tk else 0.0)
                cov_top = PF.coverage([u.allele for u in tk],
                                      universe if arm == "mhc1" else None)
                o.row("arm", d["n_chosen"], arm,
                      f"{d['p_at_least']:.4f}",
                      f"target >={d['target']} of {d['slots']} slots; "
                      f"top-{d['slots']}-by-score gives {p_top:.4f}; "
                      f"H/Hmax {d['coverage']['entropy_ratio']:.3f} vs "
                      f"{cov_top['entropy_ratio']:.3f}, gini {d['coverage']['gini']:.3f} vs "
                      f"{cov_top['gini']:.3f}")
                for u in d["units"]:
                    o.row("composed", arm, u.gene, f"{u.p:.4f}",
                          f"{u.allele} {u.kind} {u.cls}")
                for u in tk:
                    o.row("top-rank", arm, u.gene, f"{u.p:.4f}",
                          f"{u.allele} {u.kind} {u.cls}")
            o.row("composition", len(comp.units), "joint", f"{comp.joint:.4f}",
                  f"gini {comp.coverage['gini']:.3f}, H/Hmax "
                  f"{comp.coverage['entropy_ratio']:.3f} over "
                  f"{comp.coverage['n_covered']}/{comp.coverage['n_allotypes']} allotype(s)")
        for allele, (n, s, y) in sel.per_allele().items():
            o.row("allotype", n, allele, f"{y:.4f}", f"sum p = {s:.4f}")
        for t in sel.trace:
            if not t["kept"]:
                o.row("not selected", "", t["gene"], f"{t['p']:.4f}",
                      f"{t['allele']} threshold {t['threshold']:.4f}")
        for name, c in built:
            # One cassette keeps the plain section names it has always had; two qualify them, so a
            # reader can tell the composed layout from the score-only one without counting rows.
            pre = "" if len(built) == 1 else f"{name.split('_', 1)[-1]}:"
            for i, (u, (lo, hi)) in enumerate(zip(c.units, c.boundaries), 1):
                o.row(f"{pre}unit", i, u.gene, u.allele,
                      f"{lo}-{hi} p={u.p:.4f} {u.kind} {u.cls}")
            for j in c.junctions:
                o.row(f"{pre}junction", j["left"] + 1, f"{j['left'] + 1}|{j['right'] + 1}",
                      f"{j['score']:.4f}", f"{j['peptide']} at {j['offset']}")
            o.row(f"{pre}cassette", len(c.units), f"spacer={c.spacer}", f"{c.cost:.4f}",
                  f"{len(c.sequence)} aa, worst junction {c.worst_junction:.4f}")
            o.row(f"{pre}sequence", "", "", c.sequence, "")
    finally:
        o.close()

    if built and a.fasta:
        with open(a.fasta, "w") as fh:
            for name, c in built:
                fh.write(f">{name} units={len(c.units)} spacer={c.spacer} "
                         f"objective={a.objective} n0={a.n0:g}\n{c.sequence}\n")
        print(f"# wrote {a.fasta}: {len(built)} cassette(s)", file=sys.stderr)

    if cas and (a.map_tsv or a.map_json):
        alleles1 = _read_alleles(a.alleles) or sorted({u.allele for u in cas.units if u.allele})
        s1 = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=("mhc1",))
        r1 = V.store_ranker(s1, alleles1, cls="mhc1")
        r2 = None
        a2 = _read_alleles(a.map_alleles_mhc2)
        if a2:
            s2 = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=("mhc2",))
            r2 = V.store_ranker(s2, a2, cls="mhc2")
        else:
            print("# no --map-alleles-mhc2: the map is class I only, so `self_help` -- whether a "
                  "unit's CD8 epitope has overlapping CD4 help from the same unit -- is not "
                  "computed", file=sys.stderr)
        feats = V.epitope_map(cas, r1, r2, threshold=a.map_threshold)
        summ = V.write_map(cas, feats, tsv_path=a.map_tsv, json_path=a.map_json)
        print(f"# map: {summ['n_mhc1']} class-I and {summ['n_mhc2']} class-II epitope(s) over "
              f"{summ['length_aa']} aa, {summ['n_junction_spanning']} spanning a junction; "
              f"{summ['n_units_with_self_help']}/{summ['n_units']} unit(s) carry their own "
              f"class-II help", file=sys.stderr)
        for path in (a.map_tsv, a.map_json):
            if path:
                print(f"# wrote {path}", file=sys.stderr)

    if built and a.fasta_nt:
        with open(a.fasta_nt, "w") as fh:
            for name, c in built:
                cds = V.back_translate(c.sequence)
                fh.write(f">{name}_cds units={len(c.units)} spacer={c.spacer} "
                         f"nt={len(cds)}\n{cds}\n")
                print(f"# {a.fasta_nt}: {name} {len(cds)} nt, "
                      f"{len(V.slippery_sites(cds))} slippery site(s) remaining", file=sys.stderr)
        print(f"# wrote {a.fasta_nt}", file=sys.stderr)


def cmd_deslip(a):
    """Scan a cassette CDS for the m1-pseudouridine +1 frameshift motif, and optionally repair it."""
    from . import vector as V

    cds = _read_seq(a.cds)
    sites = V.slippery_sites(cds)
    o = _Out(a, "site")
    try:
        o.header("codon_index", "nt_offset", "codon", "next_codon")
        for s in sites:
            o.row(s["codon_index"], s["nt_offset"], s["codon"], s["next_codon"])
    finally:
        o.close()
    print(f"# {len(sites)} slippery site(s) in {len(cds)} nt "
          f"({len(cds) // 3} codons)", file=sys.stderr)
    if not sites:
        print("# nothing to repair. Note this motif only matters for an m1-pseudouridine construct; "
              "on unmodified uridine the scan is not the relevant check", file=sys.stderr)
    if a.fix:
        fixed, n = V.deslip(cds)
        with open(a.fix, "w") as fh:
            fh.write(f">deslipped n_fixed={n}\n{fixed}\n")
        print(f"# wrote {a.fix}: {n} codon(s) rewritten TTT -> TTC, protein unchanged",
              file=sys.stderr)


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
    _add_mhc2_report(r)

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
    _add_mhc2_report(bd)

    s = sub.add_parser("scan", help="find presented peptides in a protein (sequence or FASTA path)")
    s.add_argument("protein")
    s.add_argument("--allele")
    s.add_argument("--cls", choices=("mhc1", "mhc2"))
    s.add_argument("--top", type=int, default=3)
    s.add_argument("--correction", choices=("bonferroni", "bh"),
                   help="multiple-testing control over windows x alleles (FWER / BH-FDR)")
    _add_store_opts(s)
    s.set_defaults(fn=cmd_scan)
    _add_mhc2_report(s)

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
    pr.add_argument("--core", action="store_true",
                    help="append the NetMHCpan-style 9-residue binding core, its 0-based "
                         "offset and the register behind it. Distinct from --footprint core, "
                         "which changes what the model scores; this only reports")
    _add_store_opts(pr)
    pr.set_defaults(fn=cmd_predict)
    _add_mhc2_report(pr)

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
    rk.add_argument("--score", choices=("aggregate", "gate"), default="aggregate",
                    help="`aggregate` (default) scores with the fitted model in "
                         "data/aggregate_mhc1.json -- the one the benchmark fitted. `gate` is the "
                         "two-term noisy-AND that was the default before 0.19.0, kept so a run can "
                         "be compared against the old ordering")
    rk.add_argument("--extended", action="store_true",
                    help="append the mimicry aggregate: signed viral / self / thymus contributions "
                         "per anchor and TCR-facing channel, their sum, and the autoimmunity "
                         "read-out. Appended as COLUMNS -- the ranking is unchanged, because "
                         "whether mimicry belongs inside the gate is not settled")
    rk.add_argument("--annotate", action="store_true",
                    help="append what each candidate actually resembles: the nearest self / viral / "
                         "thymic mimic per channel with its source protein, plus the nearest "
                         "validated neoantigen and its distance")
    rk.add_argument("--core", action="store_true",
                    help="append the NetMHCpan-style 9-residue binding core, its 0-based offset "
                         "and the register behind it: for class I both anchors are held and the "
                         "central bulge gives way (an 8-mer is its own core), for class II the "
                         "register-anchored 9-mer. Reported, never scored")
    rk.add_argument("--no-self", action="store_true",
                    help="with --extended/--annotate, skip the host proteome. It is the expensive "
                         "reference (several minutes, ~7 GB) and carries the largest coefficients, "
                         "so this is faster and deliberately a smaller model")
    rk.add_argument("--prevalence", type=float, default=round(37.0 / 615.0, 4),
                    help="assumed fraction of this candidate pool that responds, used to put "
                         "`score` on a probability axis as the `p_response` column (default "
                         "%(default)s -> TESLA's 37 of 615). It is a PRIOR you own: the fit gave "
                         "every screen its own intercept precisely so base rate stayed out of the "
                         "slopes, and the nine screens behind it span 0.0060%% to 59.7%% positive. "
                         "It shifts every probability and moves no rank")
    rk.add_argument("--top", type=int, help="print only the top N candidates")
    rk.add_argument("--out", help="write TSV here instead of stdout")
    _add_store_opts(rk)
    rk.set_defaults(fn=cmd_rank)
    _add_mhc2_report(rk)

    ex = sub.add_parser("explain", help="every component of the aggregate for one (peptide, allele)")
    ex.add_argument("peptide", nargs="?", default="")
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

    ng = sub.add_parser("neoag",
                        help="annotate candidates against the tested-neoantigen database "
                             "(nearest validated-immunogenic peptide + substitution distance)")
    ng.add_argument("input", nargs="*", help="peptide(s); or use --peptides")
    ng.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ng.add_argument("--max-subs", type=int, default=2,
                    help="fuzzy search radius. 0 is an exact database lookup; 2 roughly doubles "
                         "the recall of a fresh cohort's true positives over exact lookup")
    ng.add_argument("--known-only", action="store_true",
                    help="keep only exact database matches (distance 0)")
    ng.add_argument("--hits-only", action="store_true",
                    help="drop candidates with nothing inside --max-subs")
    ng.add_argument("--core", action="store_true",
                    help="append the binding core, its 0-based offset and the register behind it. "
                         "This command has no allele, so a class-II core is the allele-agnostic "
                         "register and `core_source` reads `heuristic`")
    _add_batch_opts(ng, "candidate")
    ng.set_defaults(fn=cmd_neoag)

    my = sub.add_parser("mimicry",
                        help="the fitted mimicry aggregate: signed viral / self / thymus "
                             "contributions per anchor and TCR-facing channel, and their sum")
    my.add_argument("input", nargs="*", help="peptide(s); or use --peptides")
    my.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    my.add_argument("--corpus", nargs="?", const="screens", default=None,
                    help="map the aggregate to a probability against a NAMED fitted corpus "
                         "(default 'screens'). Omit to keep the prior-free log-odds, which is what "
                         "you rank on -- the screens run from 0.048%% to 46.8%% positive, so an "
                         "unqualified probability mostly reports one of those prevalences")
    my.add_argument("--annotate", action="store_true",
                    help="also emit what was hit: the nearest mimic per channel and its source "
                         "protein. Without this the table is the eight numeric columns")
    my.add_argument("--no-self", action="store_true",
                    help="skip the host proteome. It is the expensive reference and it carries the "
                         "largest coefficients, so this scores a deliberately smaller model")
    my.add_argument("--coefficients", action="store_true",
                    help="print the shipped model's coefficients and fit record; score nothing")
    _add_batch_opts(my, "candidate")
    my.set_defaults(fn=cmd_mimicry)

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

    vc = sub.add_parser("vector",
                        help="assemble a polyepitope cassette: withdraw on safety, choose how many "
                             "units per allotype, order them, pick a spacer")
    vc.add_argument("--candidates", required=True, metavar="FILE",
                    help="TSV of units: peptide, gene, allele, p (+ optional mutation_index, cls). "
                         "`peptide` is the LONG window around the mutation, not the minimal epitope "
                         "-- a minimal peptide loads onto any cell without costimulation and is the "
                         "tolerising configuration. `-` = stdin")
    vc.add_argument("--context", metavar="FILE",
                    help="the window FASTA `rank` was run on. With it, --candidates may be `rank`'s "
                         "own output of MINIMAL epitopes: each is joined back to its source window "
                         "and re-centred as a long unit, one per variant rather than one per "
                         "register. Without it --candidates must already carry long windows")
    vc.add_argument("--unit-length", type=int, default=27, metavar="N",
                    help="unit window length for --context (default 27, the BioNTech backbone "
                         "configuration; see mhcmatch.vector.unit)")
    vc.add_argument("--n0", type=float, required=True, metavar="F",
                    help="per-allotype capacity, the one free parameter of the stopping rule. "
                         "REQUIRED and with no default on purpose: nothing in the public record fits "
                         "it, so the value is yours to defend and it is recorded in the output")
    vc.add_argument("--quota", metavar="ARM=SLOTS:TARGET",
                    help="compose the cassette to quotas instead of taking the ranked top, e.g. "
                         "'mhc1=8:2,mhc2=4:1,nonconventional=3:1' -- eight class-I slots of which "
                         "at least two should respond, and so on. Arms are disjoint: a unit whose "
                         "`kind` column is anything but `missense` is charged to `nonconventional`, "
                         "so the constraint bites. The same slot budgets filled by score alone are "
                         "reported beside it, because 'not the same as ranking' is a claim about "
                         "YOUR candidates")
    vc.add_argument("--block-live", type=float, default=0.5, metavar="Q",
                    help="P(a block is live) in the response model behind --quota (default "
                         "%(default)s). A block is an allotype: if the recipient never mounts a "
                         "response on that allotype, none of its units respond however good they "
                         "are. Measure it on your own readout with "
                         "mhcmatch.portfolio.betabinom_rho before trusting a default")
    vc.add_argument("--evenness", type=float, default=0.0, metavar="W",
                    help="weight on class-I allotype evenness (H/Hmax) in --quota's objective "
                         "(default %(default)s = off). The block model already prefers spread when "
                         "spread helps; this is for when it does not and you want it anyway. "
                         "Homozygosity is handled: the denominator is the DISTINCT allotypes in "
                         "--alleles, so a homozygous locus is not scored as a design flaw")
    vc.add_argument("--alleles", help="the recipient's allotypes for junction scoring "
                                      "(comma-separated or a file); default = those in the table")
    vc.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    vc.add_argument("--cls-filter", action="store_true",
                    help="select only units whose own `cls` matches --cls")
    vc.add_argument("--screen", action="store_true",
                    help="withdraw units on essential-tissue risk BEFORE selecting. Costs a "
                         "whole-proteome index (minutes, several GB); without it no safety check "
                         "runs at all and the cassette carries whatever it was handed")
    vc.add_argument("--min-tpm", type=float, default=0.25, metavar="F",
                    help="essential-tissue expression floor for --screen. 0.25 because MAGE-A12 sits "
                         "at 0.33 TPM in brain and killed two patients; a conventional 5 would pass "
                         "it")
    vc.add_argument("--max-subs", type=int, default=0, metavar="N",
                    help="self-origin search radius for --screen. 0 = exact coincidence, which is "
                         "the default because the decision is per unit while the search is per "
                         "register: at radius 1 over 8-11mers, 3 of 6 random 27-mers get withdrawn "
                         "by chance. Raise it only together with dropping 8-mers")
    vc.add_argument("--objective", default="sum", choices=("sum", "rate"),
                    help="junction cost: `sum` of the strongest binder per junction (pVACvector's "
                         "logic, biased toward the shortest spacer), or `rate` = binders per "
                         "register, which is length-neutral and needs --binder-threshold. The two "
                         "disagree on real payloads, so choose")
    vc.add_argument("--binder-threshold", type=float, metavar="F",
                    help="-log10(%%rank) above which a junction window counts as a binder; "
                         "required by --objective rate")
    vc.add_argument("--threshold", type=float, metavar="F",
                    help="stop at the first spacer whose worst junction falls at or below this, "
                         "instead of trying them all and taking the cheapest")
    vc.add_argument("--fasta", metavar="FILE", help="also write the cassette sequence as FASTA")
    vc.add_argument("--fasta-nt", metavar="FILE",
                    help="also write the cassette CODING SEQUENCE as FASTA -- highest-usage human "
                         "codon per residue, backed off to avoid homopolymers, then deslipped. "
                         "Epitope cassette only: no start, no stop, no leader, no trafficking "
                         "domain")
    vc.add_argument("--map", metavar="FILE", dest="map_tsv",
                    help="also write the cassette MAP as TSV: one row per unit, linker and "
                         "predicted epitope, with 1-based coordinates over the cassette, the "
                         "presenting allele, and which class-I and class-II epitopes overlap each "
                         "other. A peptide presented by two of the recipient's alleles gets TWO "
                         "rows -- at a heterozygous locus those are two presentation events")
    vc.add_argument("--map-json", metavar="FILE",
                    help="the same map as JSON, plus the per-unit summary and the sequence, which "
                         "is what a viewer needs to draw the cassette without recomputing anything")
    vc.add_argument("--map-threshold", type=float, default=2.0, metavar="F",
                    help="%%rank at or below which a window enters the map (default 2.0)")
    vc.add_argument("--map-alleles-mhc2", metavar="LIST",
                    help="the recipient's class-II allotypes (comma-separated or a file). Without "
                         "them the map carries class I only, and a unit's `self_help` column -- "
                         "whether its CD8 epitope has overlapping CD4 help from the SAME unit -- "
                         "cannot be computed")
    _add_store_opts(vc)
    vc.add_argument("--out", metavar="FILE", help="write the report TSV here instead of stdout")
    vc.set_defaults(fn=cmd_vector)

    ds = sub.add_parser("deslip",
                        help="find (and repair) the m1-pseudouridine +1 frameshift motif in a "
                             "cassette coding sequence")
    ds.add_argument("cds", help="coding sequence, or a FASTA path (T or U, case-insensitive)")
    ds.add_argument("--fix", metavar="FILE",
                    help="write the repaired CDS here: every TTT before a T/C-starting codon becomes "
                         "TTC, which is synonymous, so the protein is unchanged")
    ds.add_argument("--out", metavar="FILE", help="write the site TSV here instead of stdout")
    ds.set_defaults(fn=cmd_deslip)

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
