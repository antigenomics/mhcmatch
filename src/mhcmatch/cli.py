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
``mimics``, ``genes``) releases the GIL and scales across cores. The scoring heads are small numpy
products per peptide, so threads there would buy nothing and the flag is not offered rather than
being offered and ignored.
"""
from __future__ import annotations

import argparse
import contextlib
import gzip
import os
import re
import sys
import time

from . import Proteome, Store, __version__, pseudoseq
from .cassette import RHO_ASSAYED as CA_RHO
from .rank import POOL_PREVALENCE

#: Verbosity, set once by :func:`main` from ``-v`` / ``-q``. 0 quiet, 1 normal, 2 verbose.
#: Module-level because every ``cmd_*`` reads it and threading it through 19 signatures would be
#: ceremony -- the CLI is a single process with one user.
_V = 1


def say(msg: str, level: int = 1, flush: bool = False) -> None:
    """One progress line to **stderr**, ``#``-prefixed.

    Stderr, always: stdout is a TSV or FASTA stream that callers pipe, and a progress line in it is
    a corrupt row. The ``#`` prefix is the house convention and predates this helper -- what the
    helper adds is that the stream and the prefix can no longer disagree, which they did in 25 of
    51 places.

    ``level=2`` is ``-v`` only. Use it for wall clock and for anything a normal run should not
    narrate; ``level=1`` is what a user watching a multi-minute command deserves to see.
    """
    if _V >= level:
        print(f"# {msg}", file=sys.stderr, flush=flush)


@contextlib.contextmanager
def step(label: str, level: int = 1):
    """Announce an expensive step before it runs, and time it after.

    The announcement is normal-verbosity and the elapsed line is ``-v``: a user staring at a silent
    terminal needs to know *what* is taking the time without being asked to opt in, but the
    duration is diagnostics. Nine of nineteen subcommands used to print nothing at all through
    multi-minute work -- ``binder`` sat through a ~45 s calibrator build in silence with the cost
    recorded only in a source comment.

    ``level=2`` makes the whole thing ``-v`` only, which is what the outer wrapper in :func:`main`
    uses -- announcing every invocation of an instant command like ``decompose`` is noise.
    """
    say(f"{label} ...", level, flush=True)
    t = time.perf_counter()
    yield
    say(f"{label}: {time.perf_counter() - t:.1f} s", max(level, 2))


def _add_verbosity(p) -> None:
    """``-v`` / ``-q`` on a subparser, so both ``mhcmatch -v vector`` and ``mhcmatch vector -v``
    work. Argparse binds a top-level flag only *before* the subcommand, and the second form is the
    one people type."""
    g = p.add_mutually_exclusive_group()
    g.add_argument("-v", "--verbose", action="store_true",
                   help="per-step wall clock and extra detail, on stderr")
    g.add_argument("-q", "--quiet", action="store_true",
                   help="suppress progress; errors and real output are unaffected")


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
    with step(f"loading the {a.species} pmhc panel ({a.tier})"):
        return Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species)


@contextlib.contextmanager
def _open_text(path):
    """Read ``path`` as text --- gzipped if it ends ``.gz``, stdin if it is ``-``.

    Four readers below opened their input with the same five lines and the same ``try``/``finally``
    that must not close stdin. Only the opening is shared: the parsing and the error messages stay
    in each reader, because those are the part that differs and the part a caller reads when their
    table is wrong.
    """
    fh = sys.stdin if path == "-" else (
        gzip.open if str(path).endswith(".gz") else open)(path, "rt")
    try:
        yield fh
    finally:
        if path != "-":
            fh.close()


def _read_seq(arg):
    """A raw sequence, or the concatenated sequences of a FASTA file path."""
    if os.path.exists(arg):
        from .proteome import read_fasta
        seqs = read_fasta(arg)
        if seqs:
            return "".join(seqs.values())
    return arg.strip()


#: How a table may spell its peptide column, best first. **One constant, because every reader in
#: this module needs the same answer**: `_read_peptides`, `_read_table` and `_cassette_rows` each
#: resolved it separately, and a pipeline candidate table -- which spells it `epitope` -- was
#: therefore accepted by `rank` and refused by `neoag`, `mimicry` and `cassette select` in the same
#: chain. A caller should not have to rename a column between two of our own commands.
PEPTIDE_COLUMNS: tuple = ("peptide", "epitope")
#: The per-unit response probability, under either spelling. `rank` writes `p_response`;
#: `cassette build/order` read `p`. They were two names for one number, and the error message
#: for the mismatch told the caller to rename the column by hand -- which made the worked
#: example in the README exit 1 as written. Resolved like :data:`PEPTIDE_COLUMNS`: the first
#: spelling present wins, and `p` is what the row carries afterwards.
RESPONSE_COLUMNS: tuple = ("p", "p_response")


def _read_peptides(path, inline=()):
    """Peptides from ``path`` (one per line, or the ``peptide`` column of a TSV) plus any inline.

    ``-`` reads stdin, so this composes with a pipe. The column may also be spelled ``epitope``
    (:data:`PEPTIDE_COLUMNS`). Whole-file reads on purpose: the scoring paths are vectorised or
    amortised over one setup, so handing them a whole deposit is both the fast path and the
    intended one."""
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
        name = next((c for c in PEPTIDE_COLUMNS if c in cols), None)
        col = cols.index(name) if name else None
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
    with _open_text(path) as fh:
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
    if getattr(a, "out", None) or getattr(a, "tsv", False):
        # One row per (window, allele): the aligned form collapses the alleles into one cell, which
        # cannot be grouped or counted downstream.
        say(f"{len(hits)} presented window(s){label}")
        out = _Out(a, "hit")
        out.header("position", "peptide", "allele", "enrichment", "n_votes", "binder")
        for pos, pep, binders in hits:
            for b in binders:
                out.row(pos, pep, _allele(a, b.allele), f"{b.enrichment:.6g}", b.n_votes,
                        int(bool(b.binder)))
        out.close()
        return
    print(f"# {len(hits)} presented window(s){label}")
    for pos, pep, binders in hits:
        print(f"{pos:>5}  {pep:<14}  {','.join(_allele(a, b.allele) for b in binders)}")


def cmd_source(a):
    with step(f"loading the {a.proteome} proteome"):
        pm = (Proteome.from_fasta(a.proteome) if os.path.exists(a.proteome)
              else Proteome.from_hf(a.proteome))
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
        say(f"wrote {a.native}")
    if a.scored_csv:
        P.write_scored_csv(preds, a.scored_csv, core=a.core)
        say(f"wrote {a.scored_csv}")
    if not a.native and not a.scored_csv:
        print(f"# {len(preds)} predicted binder(s) (%rank <= {a.rank_threshold}) over "
              f"{len(alleles)} allele(s)")
        for p in preds[:(a.top or 20)]:
            print(f"{p.peptide:<15} {_allele(a, p.allele):<18} %rank={p.percent_rank:<6} {p.band:<11} "
                  f"{p.var.get('gene_name', '')}")


def cmd_logo(a):
    from . import logo
    m = logo.motif(_store(a), a.allele, a.cls or "mhc1")
    if getattr(a, "out", None) or getattr(a, "tsv", False):
        # The full PWM, one row per (position, residue) -- the aligned form keeps only the top 3,
        # which is a display choice and not something a figure should inherit.
        say(f"{a.allele} width={m['width']} n={m['n']} "
            f"lengths={dict(sorted(m['length_hist'].items()))}")
        out = _Out(a, "cell")
        out.header("pos", "bits", "aa", "p")
        for i, (bits, col) in enumerate(zip(m["bits"], m["pwm"]), 1):
            for aa, q in sorted(col.items(), key=lambda x: -x[1]):
                out.row(i, f"{bits:.4f}", aa, f"{q:.6g}")
        out.close()
        return
    print(f"# {a.allele}  width={m['width']}  n={m['n']}  lengths={dict(sorted(m['length_hist'].items()))}")
    for i, (bits, col) in enumerate(zip(m["bits"], m["pwm"]), 1):
        top = sorted(col.items(), key=lambda x: -x[1])[:3]
        print(f"  pos {i:>2}  {bits:4.2f} bits  " + " ".join(f"{aa}:{p:.2f}" for aa, p in top))


def _read_alleles(arg):
    """``--alleles`` accepts a comma-separated list or a path to a file holding one (pipeline form)."""
    if arg and os.path.exists(arg):
        arg = open(arg).read().strip()
    return [x.strip() for x in (arg or "").replace("\n", ",").split(",") if x.strip()]


#: Loci a class-I panel is built from, and the class-II beta loci a pair key is keyed on. A typing
#: file lists both classes in one table, and the locus is the only thing that says which is which.
_MHC1_LOCI = ("A", "B", "C")
_MHC2_BETA = ("DRB1", "DRB3", "DRB4", "DRB5", "DPB1", "DQB1")
_MHC2_ALPHA = {"DPB1": "DPA1", "DQB1": "DQA1"}


def _typing_rows(path):
    """Allele names out of a typing file, whatever shape it is in.

    Three shapes, because three tools write them: a TSV with an ``Allele`` (or ``allele``) column --
    what OptiType, kourami and HLA-LA emit and what a donor's own ``.alleles.tsv`` is; a
    comma-separated list; one name per line. A ``Locus`` column is not required: the locus is read
    off the name, which is the only thing present in all three.
    """
    text = open(path).read() if os.path.exists(path) else (path or "")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if lines and "\t" in lines[0]:
        cols = lines[0].rstrip("\n").split("\t")
        col = next((c for c in cols if c.strip().lower() == "allele"), None)
        if col is not None:
            i = cols.index(col)
            return [f.split("\t")[i].strip() for f in lines[1:] if len(f.split("\t")) > i]
    return [x.strip() for ln in lines for x in ln.split(",") if x.strip()]


def _locus_of(name):
    """``'HLA-DQB1*03:01'`` -> ``'DQB1'``, ``'A*01:01'`` -> ``'A'``; ``''`` when it is not an HLA name."""
    n = re.sub(r"^HLA[- ]", "", (name or "").strip(), flags=re.I)
    m = re.match(r"^(D[PQR][AB]\d|[ABC])w?[*\d]", n, flags=re.I)
    return m.group(1).upper() if m else ""


def cmd_alleles(a):
    """A donor's typing file -> the allele list every other command's ``--alleles`` accepts.

    Three things stand between a typing file and a scored run, and each of them fails **silently**:

    * **Field depth.** Every typer writes three or four fields and a G-group suffix
      (``A*01:01:01G``); the pseudosequence tables are keyed at two. See
      :func:`mhcmatch.pseudoseq.trim_allele`.
    * **The class split.** One file lists both classes, and a class-I panel handed a DQB1 name
      resolves it to nothing.
    * **The DP/DQ pairing.** A DP or DQ molecule is an alpha-beta heterodimer and its key names both
      chains, so the two rows of a typing file have to be *joined* -- ``DQA1*05:01`` alone is not a
      molecule and does not resolve. :func:`mhcmatch.pseudoseq.class2_key` is the join; where the
      alpha is absent it is imputed from the beta (:func:`mhcmatch.pseudoseq.alpha_prior`), which is
      what DR needs too since DRA is monomorphic.

    Every dropped name is reported, because :meth:`mhcmatch.store.Store._allele_set` drops what it
    cannot find without saying so -- the whole failure mode this command exists to make loud.
    """
    from .pseudoseq import class2_key, resolve_allele, trim_allele

    names = [trim_allele(n) for n in _typing_rows(a.input)]
    by_locus = {}
    for n in names:
        loc = _locus_of(n)
        if loc and n not in by_locus.setdefault(loc, []):
            by_locus[loc].append(n)

    out, dropped = [], []
    if a.cls == "mhc1":
        for loc in _MHC1_LOCI:
            out += by_locus.get(loc, [])
    else:
        for beta_loc in _MHC2_BETA:
            alpha_loc = _MHC2_ALPHA.get(beta_loc)
            alphas = by_locus.get(alpha_loc, []) if alpha_loc else []
            for beta in by_locus.get(beta_loc, []):
                # No alpha typed -> one key with the imputed alpha. Two alphas typed -> both pairings,
                # because a heterozygous DQA1 with a heterozygous DQB1 really can present as four
                # molecules and the typing does not say which trans pairs form.
                out += [class2_key(al, beta) for al in alphas] or [class2_key("", beta)]
    seen, keep = set(), []
    for n in out:
        key, exact = resolve_allele(n, a.cls)
        if key is None:
            dropped.append(n)
        elif key not in seen:
            seen.add(key)
            keep.append(n if a.form == "input" else key)
    unknown = [n for n in names if not _locus_of(n)]
    for n in unknown:
        dropped.append(n)
    if dropped:
        say(f"dropped {len(dropped)} name(s) that resolve to no pseudosequence: "
            + ", ".join(sorted(set(dropped))), level=1)
    # A mouse haplotype is a property of the inbred line, so there is no typing file for it and this
    # command has no locus grammar that would match one. Say that, rather than leaving a run to
    # discover an empty allele list on its own -- which is exactly the silence this command exists
    # to break.
    if not keep and any(n.upper().startswith(("H-2", "H2-", "I-")) for n in names):
        say("every name looks like a mouse H-2 allele: this command reads HLA typing files. An "
            "inbred line's haplotype is a property of the line, so pass it directly, e.g. "
            "--alleles 'H2-K*d,H2-D*d,H2-L*d' / --alleles_mhc2 'H-2-IAd,H-2-IEd'", level=1)
    say(f"{len(keep)} {a.cls} allele(s) from {len(names)} typed name(s)", level=1)
    text = ",".join(keep)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text + "\n")
    else:
        print(text)


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

    **There is no search here at all.** The channels are a
    :func:`mhcmatch.mimicry.corpus_spectrum` table contraction rather than a neighbour search, so
    all three together cost three tables of 64 KB and the ranking path builds no trie -- against the
    host-proteome reference index a search would force, 6 min 15 s and ~7.5 GB.
    ``--no-self`` still matters for ``--extended``/``--annotate``,
    which report *which* reference peptide was nearest and do need the index, and for the safety
    scan.

    ``species`` keys every corpus channel, not only ``self``: all six ``thymus``/``self``/``viral``
    tables ship for each of mouse and human, so a mouse run is scored against mouse references
    throughout. The mouse thymic and viral tables are *thinner* than their human counterparts --
    25,264 and 40,244 reference windows against 140,482 and 136,618 for class I -- which is a
    precision statement about a mouse score, not a substitution. This used to read "the thymus and
    viral deposits are human-only", which stopped being true when those tables shipped.

    The ``C_phys`` pair is deliberately absent: :func:`mhcmatch.rank._finish` computes both, because
    they are matrix products against published residue vectors and need no deposit at all.
    """
    from . import mimicry as MM
    from . import rank as R

    def channels(peptides):
        # `k`, the face mask and the substitution kernel all come from the artifact
        # (`MM.corpus_geometry`), exactly as `kappa` already did. They are three halves of one
        # definition: a `kappa` fitted against a graded kernel scored under the Hamming one is a
        # different feature, not a smaller effect.
        g = MM.corpus_geometry()
        spec = MM.corpus_spectrum(cls=cls, components=("thymus", "self", "viral"),
                                  k=g["k"], self_species=species, mask=g["mask"],
                                  kernel=g["kernel"])
        rows = MM.corpus_R(list(peptides), spec, cls=cls)
        return {f"C_corpus_{c}": [r.get(c, float("nan")) for r in rows]
                for c in ("thymus", "self", "viral")
                if f"C_corpus_{c}" in R.CHANNEL_COLUMNS}

    return channels


def _rank_model(a):
    """Print the shipped aggregate itself instead of scoring anything.

    Both tables are read out of ``data/aggregate_mhc1.json`` -- the artifact the benchmark fitted
    and the library ships. Nothing is refitted here, so a figure built on this output and a run of
    ``rank`` are the same model by construction rather than by a comparison someone has to make."""
    from . import rank as R
    m = R.aggregate()
    f = m["fit"]
    say(f"{m['model']} v{m['version']}, fitted by {m['generator']}")
    say(f"{f['rows']:,} rows / {f['positives']:,} positives over {len(f['screens'])} screens; "
        f"BIC {f['bic']:.1f}, ridge tau {f['tau']}")
    say(f"intervals from {f['n_boot']} resamples of {f['bootstrap_unit']}; holdout {f['holdout']}")
    say("no global intercept: every screen was given its own, unpenalised")
    out = _Out(a, "term" if not getattr(a, "holdout", False) else "screen")
    if getattr(a, "holdout", False):
        out.header("screen", "n", "pos", "neg", "auroc", "decided")
        for r in m["loo"]:
            out.row(r["level"], r["n"], r["pos"], r["neg"], f"{r['auroc']:.4f}",
                    "yes" if r["decided"] else "no")
        for k in ("cv_peptide", "cv_twin"):
            c = m[k]
            out.row(k, "", "", "", f"{c['median_decided']:.4f}", f"{c['folds']}-fold, "
                    f"{c['groups']:,} groups, {c['n_decided']} decided, "
                    f"pooled {c['pooled_auroc']:.4f}")
        v = m["verdict"]
        say(f"verdict: {v['improvements']} improvement(s), {v['ties']} tie(s), "
            f"{v['regressions']} regression(s)")
    else:
        block = {t: b for b, ts in m["blocks"] for t in ts}
        out.header("block", "term", "coef", "sd", "boot_sd", "z", "p",
                   "ci_low", "ci_high", "sign_stab")
        for i, t in enumerate(m["features"]):
            lo, hi = m["ci95"][i]
            out.row(block.get(t, ""), t, f"{m['coef'][i]:+.4f}", f"{m['sd'][i]:.4f}",
                    f"{m['boot_sd'][i]:.4f}", f"{m['z'][i]:+.2f}", f"{m['p'][i]:.3g}",
                    f"{lo:+.4f}", f"{hi:+.4f}", f"{m['sign_stability'][i]:.4f}")
        say("blocks are entered in pipeline order, each on top of the last: a recognition "
            "coefficient is what it is worth AFTER presentation and expression")
    out.close()


def cmd_rank(a):
    """Rank neoantigen candidates from a window FASTA or an already-scored table.

    With ``--score aggregate`` (the default) every one of the model's features is computed
    *before* scoring and emitted as a column -- see :func:`_aggregate_channels`. A model emits the
    features it used and refuses to run without them.

    ``--no-self`` and ``--score aggregate`` compose: ``EPIC`` does not score on ``self_tcr``, so the
    host-proteome index is off the ranking path entirely.
    """
    from . import rank as R
    if getattr(a, "coefficients", False) or getattr(a, "holdout", False):
        _rank_model(a)
        return
    if not a.mode or not a.input:
        raise SystemExit("rank needs a mode and an input, or --coefficients / --holdout")
    # None -> mhcmatch.known's built-in sets; --no-known-refs -> {} -> lookup off
    refs = _load_refs(getattr(a, "refs", None)) if getattr(a, "refs", None) else \
        ({} if getattr(a, "no_known_refs", False) else None)
    carry = []            # the caller's own columns, in the caller's own order
    if a.mode == "pairs":
        store = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
        recs = _read_table(a.input)
        # The caller's OWN columns, read off the header rather than off a row: `_read_table` adds a
        # `peptide` key when the header spelled it `epitope`, and --passthrough emits the caller's
        # table, not one this command widened by a column they never sent.
        with _open_text(a.input) as _fh:
            carry = _fh.readline().rstrip("\n").split("\t")
        # **Two required columns, and everything else is the caller's business.** The peptide is
        # checked by `_read_table` above; the allele is checked here rather than discovered as an
        # empty `allele` field several minutes into a scoring run, where it reads as "this candidate
        # named no allele we know" -- a real and different state that `_unscored` handles, and one
        # a missing COLUMN should not be mistaken for.
        if recs and not any(c in carry for c in ("allele", "best_allele")):
            raise SystemExit(f"{a.input}: no `allele` / `best_allele` column (found {carry}). "
                             "`peptide`/`epitope` and one of these two are the only columns this "
                             "command requires; `wt_peptide`, `gene`/`gene_name` and `tpm` are used "
                             "when present, and every other column is carried through untouched")
        # `--context` supplies the germline arm the table does not carry. Before ranking, because
        # the wild type is what `binder_ranks` is asked for alongside the mutant.
        if getattr(a, "context", None):
            n = R.wt_from_windows(recs, a.context)
            say(f"--context: recovered a wild type for {n:,} of {len(recs):,} row(s) from "
                f"{a.context}", level=1)
        rows = R.rank_pairs(store, recs, cls=a.cls,
                            tissue=a.tissue, tumor=a.tumor, refs=refs, score=a.score,
                            prevalence=a.prevalence,
                            channels=_aggregate_channels(a.cls, a.no_self, a.species)
                            if a.score == "aggregate" else None)
    elif a.mode == "fasta":
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
    # `rank` floats an exact known-epitope match to the top of its *listing* -- a display choice,
    # documented on `Ranked.rank`, which the `rank` column does not follow. Under --passthrough the
    # file IS the caller's table re-ordered by our score, so the listing order and the rank have to
    # agree or the deliverable is sorted by something nobody asked for. Scoped to the flag: every
    # other caller's output is byte-identical.
    if getattr(a, "passthrough", False):
        rows = sorted(rows, key=lambda r: r.rank)
    rows = rows[:a.top] if a.top else rows
    cols = list(R.BASE_COLUMNS)
    if a.score == "aggregate":
        cols += list(R.EXPR_COLUMNS) + list(R.AGGREGATE_COLUMNS)
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
    # `--passthrough`: the caller's table comes back annotated and re-ordered, not replaced. Its
    # columns lead, in its own order; ours follow under `--prefix`. A join would not reproduce this
    # -- `rank` splits a multi-allele cell and the best presenter stands for the row, so the output
    # shares neither length nor allele column with the input.
    head = cols
    if getattr(a, "passthrough", False):
        pre = getattr(a, "prefix", "") or ""
        # **A name collision is an error, not a warning.** The contract is that the caller's columns
        # come back untouched and ours are appended, and two columns under one name breaks it in the
        # worst way available: every reader that keys a row by name -- `csv.DictReader`, pandas,
        # polars, our own `_read_table` -- silently resolves the duplicate in favour of one of them,
        # and which one is not something the file records. So the run stops here instead.
        clash = sorted(set(carry) & {pre + c for c in cols})
        if clash:
            raise SystemExit(
                f"--passthrough: {len(clash)} of your column(s) collide with the names this "
                f"command adds under --prefix {pre!r}: {', '.join(clash)}. Your columns are "
                "emitted unchanged and ours are appended, so the two cannot share a name -- "
                "choose a --prefix that does not collide (the shipped deliverables use `mm_`), "
                "or rename those columns upstream")
        head = carry + [pre + c for c in cols]
    out = open(a.out, "w") if a.out else sys.stdout
    try:
        print("\t".join(head), file=out)
        for i, r in enumerate(rows, 1):
            cells = [str(r.rank), r.peptide, _allele(a, r.allele),
                     _allele(a, r.allele_scored) if r.allele_scored else "",
                     r.gene, f"{r.score:.6g}",
                     f"{r.p_response:.4g}",
                     f"{r.presentation:.4g}", f"{r.binder:.4g}", f"{r.occupancy:.4g}",
                     f"{r.d_occupancy:.4g}", "1" if r.wt_absent else "0",
                     f"{r.agretopicity:.4g}",
                     f"{r.physchem:.4g}", f"{r.expression:.4g}", f"{r.expr_pct:.4g}",
                     "1" if r.expression_imputed else "0",
                     str(r.n_alleles_presenting), r.alleles_presenting,
                     r.imputed, r.wt_peptide,
                     r.known_epitope, r.variant_type]
            if a.score == "aggregate":
                # `.get`: an artifact that does not declare an expression term never sets it, and
                # an absent term is an empty cell rather than a KeyError or a fabricated 0.
                cells += ["" if r.components.get(c) is None else f"{r.components[c]:.6g}"
                          for c in R.EXPR_COLUMNS]
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
            if getattr(a, "passthrough", False):
                cells = [str(r.row.get(c, "")) for c in carry] + cells
            print("\t".join(cells), file=out)
    finally:
        if a.out:
            out.close()
            say(f"wrote {a.out}: {len(rows):,} candidate(s)")
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


def _read_table(path, col=PEPTIDE_COLUMNS):
    """Every row of a TSV with a ``peptide`` column, as dicts, preserving column order.

    ``col`` names that column for the one caller whose table does not spell it ``peptide``
    (``genes``, which offers ``--peptide-col``); it is the column normalised to upper case, so the
    key a lookup is built on and the cell that is written back agree. It may be a tuple of accepted
    spellings (default :data:`PEPTIDE_COLUMNS`), in which case the first the header carries wins --
    which is how a pipeline table spelling it ``epitope`` is read without a rename stage.

    **A row resolved to a spelling in :data:`PEPTIDE_COLUMNS` always gets a ``peptide`` key**, so a
    downstream reader never has to know which of them arrived; the caller's own column is untouched
    and still in the dict, so a passthrough that carries "every column but ``peptide``" still
    carries it. A caller-named column (``--peptide-col mt_peptide``) gets no such alias --
    ``genes`` writes its header from the row's keys and an invented one would appear in the
    output."""
    with _open_text(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        want = (col,) if isinstance(col, str) else tuple(col)
        col = next((c for c in want if c in cols), want[0])
        if col not in cols:
            raise SystemExit(f"{path}: no `{'` / `'.join(want)}` column (found {cols})")
        out = []
        for line in fh:
            line = line.rstrip("\n")
            if line:
                d = dict(zip(cols, line.split("\t")))
                d[col] = (d.get(col) or "").strip().upper()
                # Only for a spelling WE know. A caller who named their own column
                # (`genes --peptide-col mt_peptide`) gets no invented key -- `genes` writes its
                # header from the row's keys, so one would appear as a third column in the output.
                if col in PEPTIDE_COLUMNS:
                    d.setdefault("peptide", d[col])
                out.append(d)
        return out


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


def cmd_genes(a):
    """Add a ``gene`` column to a candidate table, naming the gene each peptide derives from.

    The step the expression axis cannot do without: ``expr_lvl`` and ``expr_norm`` are keyed on an
    HGNC symbol, and a deposit that ships peptides without one leaves both terms at a single
    mean-imputed constant for every row. :meth:`mhcmatch.proteome.Proteome.assign_genes` is the
    whole computation; this reads the table, carries every original column through, and writes the
    annotation back so ``rank pairs`` can be handed the result with no join.

    **A tie becomes several rows, not a refusal.** Which of several equally-near parents a peptide
    should be scored under is a question the expression reference answers and this pass cannot, so
    every tied gene is emitted and the caller takes the best aggregate score per peptide
    (``group_by(peptide).agg(max(score))``). **An unresolved peptide keeps its row** with an empty
    ``gene``, and is scored on the terms it does have -- losing the row would be the larger error.
    """
    with step(f"loading the {a.species} proteome"):
        pm = (Proteome.from_fasta(a.species) if os.path.exists(a.species)
              else Proteome.from_hf(a.species))
    rows = _read_table(a.input, a.peptide_col)
    with step(f"resolving parent genes within {a.max_subs} substitution(s)"):
        genes = pm.assign_genes([r[a.peptide_col] for r in rows], max_subs=a.max_subs,
                                threads=a.threads)
    cols = list(rows[0]) if rows else [a.peptide_col]
    if "gene" not in cols:
        cols.append("gene")
    out = _Out(a, "row")
    out.header(*cols)
    n_res = n_tie = 0
    for r in rows:
        got = genes.get(r[a.peptide_col]) or [""]
        n_res += 1 if got[0] else 0
        n_tie += 1 if len(got) > 1 else 0
        for g in got:
            out.row(*(g if c == "gene" else r.get(c, "") for c in cols))
    out.close()
    say(f"{n_res:,} of {len(rows):,} peptide row(s) resolved to a gene, {n_tie:,} with a tie; "
        "an unresolved row keeps its place with an empty `gene`")


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
    if getattr(a, "out", None) or getattr(a, "tsv", False):
        # Numeric cells, no "median "/"IQR " prefixes inside them: the aligned form is readable and
        # is not a table anything can parse.
        out = _Out(a, "row")
        out.header("key", "context", "source", "median_tpm", "q25_tpm", "q75_tpm", "n")
        rec = EX.lookup(a.key, tissue=a.tissue, tumor=a.tumor) if (a.tissue or a.tumor) else None
        if rec:
            out.row(a.key, a.tumor or a.tissue or "", rec["source"], f"{rec['median_tpm']:.6g}",
                    f"{rec['q25_tpm']:.6g}", f"{rec['q75_tpm']:.6g}", rec["n"])
        elif a.tissue or a.tumor:
            say(f"no reference row for {a.key!r} in {a.tumor or a.tissue!r}")
        if a.safety:
            for t, v in EX.safety_profile(a.key, top=a.top or 10):
                out.row(a.key, t, "GTEx", f"{v:.6g}", "", "", "")
        out.close()
        return
    rec = EX.lookup(a.key, tissue=a.tissue, tumor=a.tumor) if (a.tissue or a.tumor) else None
    if not rec:
        if a.tissue or a.tumor:
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
    # The single-pipeline GTEx/TCGA reference, as the three files scoring actually reads -- 6.6 MB
    # between them. The table they are derived from is 38.6 MB and is deliberately NOT staged here:
    # nothing on a scoring path parses it, and on a slow link it is minutes of download for numbers
    # already in the matrix. `expression.REFERENCE_TOIL_FILE` fetches it when an analysis wants rows.
    "expression/toil_matrix.npz",                 # 58,581 genes x 86 contexts, dense float32
    "expression/toil_floors.tsv",                 # per-context abundance floors, 88 rows
    "expression/context_synonyms.tsv",            # free-text origin -> context, 448 rows
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


def _cell(fields, ix, names):
    """First non-empty cell among ``names``, ``""`` if none is present. The variant class is
    spelled ``kind`` by hand, ``variant_type`` by ``rank``, and ``mm_variant_type`` by
    ``rank --prefix mm_``; ``cassette build --quota`` reads it, and reading only the first spelling
    is how a non-conventional quota comes back satisfiable by missense alone."""
    for n in names:
        if n in ix and len(fields) > ix[n] and fields[ix[n]].strip():
            return fields[ix[n]].strip()
    return ""


def _read_units(path, unit_column: str = "peptide"):
    """``[Unit]`` from a TSV with ``peptide``/``gene``/``allele``/``p`` (+ optional
    ``mutation_index``, ``cls``, ``kind``).

    ``unit_column`` names the column holding the **long window**, for a caller whose table carries
    it beside the minimal epitope rather than in place of it -- a pipeline candidate table spelling
    it ``epitope_context`` is the case in hand, and at 27 residues it is already the shipped
    ``--unit-length``. This is the alternative to ``--context``, not a second version of it:
    ``--context`` rebuilds a window from the variant's FASTA when the table has none, and this reads
    one the table already has.

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

    with _open_text(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        pcol = next((c for c in RESPONSE_COLUMNS if c in cols), None)
        need = (unit_column, "gene", "allele")
        missing = [c for c in need if c not in cols] + ([] if pcol else ["p"])
        if missing:
            raise SystemExit(f"{path}: missing column(s) {', '.join(missing)}; a unit table needs "
                             f"{', '.join(need)} and one of {' / '.join(RESPONSE_COLUMNS)} "
                             f"(+ optional mutation_index, cls). `rank` gives you gene, allele and "
                             f"a score -- `{unit_column}` must be the long window around the "
                             "mutation, not the minimal epitope")
        ix = {c: cols.index(c) for c in cols}
        ix.setdefault("p", ix[pcol])
        units = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[ix[unit_column]].strip():
                continue
            pep = f[ix[unit_column]].strip().upper()
            mi = (int(f[ix["mutation_index"]]) if "mutation_index" in ix
                  and len(f) > ix["mutation_index"] and f[ix["mutation_index"]].strip()
                  else len(pep) // 2)
            units.append(Unit(peptide=pep, mutation_index=mi, gene=f[ix["gene"]].strip(),
                              allele=f[ix["allele"]].strip(), p=float(f[ix["p"]]),
                              cls=(f[ix["cls"]].strip() if "cls" in ix and len(f) > ix["cls"]
                                   else "mhc1"),
                              kind=_cell(f, ix, ("kind", "variant_type", "mm_variant_type"))
                              or "missense"))
        return units


def _read_unit_rows(path):
    """``[dict]`` from the same TSV :func:`_read_units` reads, without the long-window contract.

    Used only with ``--context``, where ``peptide`` is deliberately a *minimal* epitope and the long
    window is rebuilt from the FASTA instead -- so the check that belongs here is that the columns
    exist, not that the peptide is long.
    """
    with _open_text(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        pcol = next((c for c in RESPONSE_COLUMNS if c in cols), None)
        need = ("peptide", "gene", "allele")
        missing = [c for c in need if c not in cols] + ([] if pcol else ["p"])
        if missing:
            raise SystemExit(f"{path}: missing column(s) {', '.join(missing)}; with --context a "
                             f"candidate table needs {', '.join(need)} and one of "
                             f"{' / '.join(RESPONSE_COLUMNS)}, which is what `rank` emits")
        rows = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(cols) or not f[cols.index("peptide")].strip():
                continue
            d = dict(zip(cols, f))
            d.setdefault("p", d[pcol])
            rows.append(d)
        return rows


def cmd_vector(a):
    """Screen, select, order: a polyepitope cassette from a table of candidate units.

    ``cassette order`` reaches the same function with ``order_only`` set, which skips the per-
    allotype sizing rule and lays out whatever units it was handed. That is the shape a caller who
    has already chosen --- with ``cassette select``, or by hand, or from a trial's published
    composition --- actually needs, and it is the same code path so the two cannot diverge.
    """
    from . import vector as V

    order_only = bool(getattr(a, "order_only", False))

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
        units = _read_units(a.candidates, getattr(a, "unit_column", None) or "peptide")
    print(f"# {len(units)} candidate unit(s) over "
          f"{len({u.allele for u in units})} allotype(s)", file=sys.stderr)

    rejected, notes, costs = [], [], {}
    if a.screen:
        from .proteome import gene_symbols
        from .store import fetch_proteome
        print(f"# screening ({a.screen_mode}): one whole-proteome window index per register length "
              f"({len(lengths)} for {a.cls}). Paid once for the whole candidate list, so screen "
              "everything in one call", file=sys.stderr, flush=True)
        fa = fetch_proteome(a.species)
        risk = V.self_origin_risk(Proteome.from_fasta(fa), gene_symbols(fa, key="accession"),
                                  min_tpm=a.min_tpm, max_subs=a.max_subs,
                                  veto_tpm=a.veto_tpm, graded=(a.screen_mode == "graded"),
                                  report_subs=a.report_subs,
                                  report_identity=a.report_identity)
        units, rejected = V.screen(units, risk, lengths=lengths, notes=notes)
        if a.report_subs and notes:
            # A d=1 coincidence is only a hazard if a T cell can see it, so the off-target's own
            # sequence has to be presented on the allotype the unit was selected for. Its own panel
            # rather than the layout binder built below: that one is scoped to the *selected* units
            # and this question is asked of every screened one, before selection has happened.
            before = len({id(u) for u, _, w in notes if "variant" in w})
            rstore = Store.from_pmhc(a.pmhc, tier=a.tier, species=a.species, classes=(a.cls,))
            ralleles = _read_alleles(a.alleles) or sorted({u.allele for u in units if u.allele})
            notes = V.presented(notes, V.store_binder(rstore, ralleles, cls=a.cls),
                                threshold=a.report_threshold)
            after = len({id(u) for u, _, w in notes if "variant" in w})
            print(f"# presentation: {before} -> {after} unit(s) carry a near-identical off-target "
                  f"that is itself presented at {10 ** -a.report_threshold:g}% rank or better",
                  file=sys.stderr)
        costs = V.offtarget_cost(notes)
        print(f"# withdrawn: {len({id(u) for u, _, _ in rejected})} unit(s), "
              f"{len(rejected)} reason(s); {len(units)} remain", file=sys.stderr)
        if notes:
            near = [n for n in notes if n[2]["clause"] == "near-identical self origin"]
            print(f"# fingerprint: {len(costs)} kept unit(s) carry {len(notes) - len(near)} "
                  f"sub-veto finding(s) below {a.veto_tpm:g} TPM, priced at --weight-offtarget "
                  f"{a.weight_offtarget:g}", file=sys.stderr)
            if near:
                print(f"# near-identity: {len({id(u) for u, _, _ in near})} kept unit(s) sit "
                      f"{a.report_subs} substitution from {len({n[2]['gene'] for n in near})} "
                      f"non-homologous expressed gene(s). Reported, not withdrawn -- read the "
                      "fingerprint before dosing", file=sys.stderr)

    if order_only:
        # Nothing to size: the caller already chose. `Selection` is still the carrier, so every
        # downstream step -- ordering, the map, the report columns -- is byte-for-byte the path a
        # selected cassette takes.
        pool = [u for u in units if not a.cls_filter or u.cls == a.cls]
        sel = V.Selection(units=pool, dropped=[], n0=float("nan"), trace=[],
                          keys=[u.allele for u in pool])
        print(f"# order only: laying out {len(sel.units)} unit(s) as given, no selection",
              file=sys.stderr)
    else:
        sel = V.select(units, n0=a.n0, cls=a.cls if a.cls_filter else None)
        print(f"# selected {len(sel.units)} of {len(units)}, expected yield "
              f"{sel.expected_yield:.2f} at n0={a.n0}", file=sys.stderr)

    comp = topk = None
    if a.quota:
        from . import portfolio as PF
        quotas = _parse_quota(a.quota)
        universe = _read_alleles(a.alleles) or sorted({u.allele for u in units if u.allele})
        comp = PF.compose(units, quotas, a.block_live, weight_evenness=a.evenness,
                          universe=universe,
                          cost=(lambda u: costs.get(u, 0.0)) if costs else None,
                          weight_cost=a.weight_offtarget)
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
                say(f"layout panel: {len(alleles)} allotype(s), paid once for "
                    f"{len(plans)} plan(s)", 2)
            with step(f"laying out {name}: {len(us)} unit(s), "
                      f"{len(us) * (len(us) - 1)} ordered pair(s) x "
                      f"{1 if a.linker else len(V.SPACERS)} spacer(s)"
                      + (f", pinned to {a.linker}" if a.linker else "")):
                built.append((name, V.order(us, binder=binder, lengths=lengths, alleles=alleles,
                                            objective=a.objective,
                                            binder_threshold=a.binder_threshold,
                                            threshold=a.threshold,
                                            linker=getattr(a, "linker", None))))
        elif us:
            # One unit has no junctions, so `order` returns before it ever calls `binder` -- building
            # the panel and its calibrators here would be ~10 s spent to lay out a cassette of one.
            built.append((name, V.order(us, binder=None, linker=getattr(a, "linker", None))))
    cas = built[0][1] if built else None

    o = _Out(a, "row")
    try:
        o.header("section", "i", "key", "value", "detail")
        # One row per (unit, register, source gene), carrying that gene's worst tissue. `screen`
        # returns a reason per tissue, and a gene is transcribed in many -- on one 7-unit test that
        # was 2,121 rows for 4 withdrawals, which is a table nobody reads. The full set is still
        # what the library returns; this is the presentation.
        def _worst(findings):
            w = {}
            for u, reg, why in findings:
                k = (id(u), reg, why.get("gene", ""))
                if k not in w or why.get("tpm", 0) > w[k][2].get("tpm", 0):
                    w[k] = (u, reg, why)
            return w.values()

        for u, reg, why in _worst(rejected):
            sub = f" {why['subs']}sub" if "subs" in why else ""
            o.row("withdrawn", "", u.gene, why.get("clause", ""),
                  f"{why.get('gene', '')}{sub} {why.get('tissue', '')} "
                  f"{why.get('tpm', 0):.1f}".strip() + (f" via {reg}" if reg else ""))
        # The off-target fingerprint of the units that were KEPT: every essential-tissue finding
        # that fell below --veto-tpm, which is the evidence a graded screen prices instead of
        # refusing. "Why was this unit kept but discounted" was previously answerable only by
        # re-running the screen.
        for u, reg, why in _worst(notes):
            o.row("fingerprint", f"{costs.get(u, 0.0):.0f}", u.gene, why.get("clause", ""),
                  f"{why.get('gene', '')} {why.get('tissue', '')} "
                  f"{why.get('tpm', 0):.2f} TPM".strip() + (f" via {reg}" if reg else ""))
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
                      f"{cov_top['gini']:.3f}"
                      + (f"; off-target cost mean {d['mean_cost']:.2f} max {d['max_cost']:.0f}"
                         if d["max_cost"] else ""))
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
                # `cassette order` does not select, so it has no --n0 and the header must not
                # claim one: formatting `None` with `:g` raised, which made `order --fasta` fail
                # for every caller who had already chosen their units.
                n0 = f" n0={a.n0:g}" if a.n0 is not None else ""
                fh.write(f">{name} units={len(c.units)} spacer={c.spacer} "
                         f"objective={a.objective}{n0}\n{c.sequence}\n")
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
        cut = V.rank_cutoffs(a.map_binder)
        t1 = a.map_threshold if a.map_threshold is not None else cut["mhc1"]
        t2 = a.map_threshold_mhc2 if a.map_threshold_mhc2 is not None else cut["mhc2"]
        mstats: dict = {}
        feats = V.epitope_map(cas, r1, r2, threshold=t1, threshold2=t2, stats=mstats)
        summ = V.write_map(cas, feats, tsv_path=a.map_tsv, json_path=a.map_json)
        print(f"# map ({a.map_binder} binders: class I %rank <= {t1:g}, class II <= {t2:g}): "
              f"{summ['n_mhc1']} class-I and {summ['n_mhc2']} class-II epitope(s) over "
              f"{summ['length_aa']} aa, {summ['n_junction_spanning']} spanning a junction; "
              f"{summ['n_units_with_self_help']}/{summ['n_units']} unit(s) carry their own "
              f"class-II help", file=sys.stderr)
        # **A zero is never left bare.** "0 class-II epitopes" cannot otherwise be told apart from
        # a ranker that never ran, a panel that resolved to nothing, and a construct whose best
        # window missed the cut by a hair. Say which, every time, and say it is a reporting cut --
        # nothing is removed from the cassette, the units table or the ranked candidates by it.
        for cls, st in sorted(mstats.items()):
            if st["n_kept"]:
                continue
            if not st["n_scored"]:
                print(f"# map: no {cls} window could be scored at all ({st['n_windows']} window(s) "
                      f"offered) -- the panel or the allele list is the thing to check, not the "
                      f"threshold", file=sys.stderr)
            else:
                print(f"# map: no {cls} epitope at %rank <= {st['threshold']:g}, and the best of "
                      f"{st['n_scored']:,} scored window(s) was %rank {st['best_rank']:.3f}. "
                      f"This is the MAP cut-off only: nothing was removed from the cassette, from "
                      f"the units table or from the ranked candidates. Widen it with "
                      f"--map-threshold{'-mhc2' if cls == 'mhc2' else ''} to annotate them",
                      file=sys.stderr)
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

    if built and a.mrna:
        # The backbone elements are read the same way alleles are: a literal, or a file holding one.
        parts = {k: _read_seq(v) if v and os.path.exists(v) else (v or "")
                 for k, v in (("leader", a.leader), ("trailer", a.trailer),
                              ("utr5", a.utr5), ("utr3", a.utr3))}
        with open(a.mrna, "w") as fh:
            for name, c in built:
                m = V.mrna(c, leader=parts["leader"], trailer=parts["trailer"],
                           utr5=parts["utr5"], utr3=parts["utr3"], poly_a=a.poly_a)
                ck = m.checks
                fh.write(f">{name}_mrna units={ck['n_units']} linker={m.linker} "
                         f"nt={ck['length_nt']} cds={ck['cds_nt']} aa={ck['protein_aa']} "
                         f"gc={ck['gc']:.3f} max_run={ck['longest_homopolymer']} "
                         f"slippery={ck['slippery_sites']} translates={ck['translates']}\n"
                         f";parts=" + ",".join(f"{q['kind']}:{q['name']}:"
                                               f"{q['start'] + 1}-{q['end']}" for q in m.parts)
                         + f"\n{m.sequence}\n")
                print(f"# {a.mrna}: {name} {ck['length_nt']} nt over {len(m.parts)} part(s), "
                      f"linker={m.linker}, GC {ck['gc']:.1%}, longest run "
                      f"{ck['longest_homopolymer']}", file=sys.stderr)
                if not ck["translates"]:
                    print(f"# WARNING {name}: the coding sequence does not read back as the "
                          f"assembled protein", file=sys.stderr)
        print(f"# wrote {a.mrna}", file=sys.stderr)


def cmd_linkers(a):
    """Print the named linker presets, so a `--linker` name never has to be guessed."""
    from . import vector as V

    o = _Out(a, "linker")
    try:
        o.header("name", "sequence", "length", "family", "class", "note")
        for name, L in V.LINKERS.items():
            o.row(name, L.sequence or "-", len(L), L.family, L.cls, L.note)
    finally:
        o.close()
    print(f"# {len(V.LINKERS)} preset(s). `class` is the class each is INTENDED for, which is "
          f"provenance and not a measurement -- pass one to --linker to pin it, or omit --linker "
          f"and let the junction scan choose", file=sys.stderr)


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


def cmd_build(a):
    """Regenerate the shipped artifacts under ``mhcmatch.data``, or just report what is stale.

    Everything mhcmatch ships is rebuilt from here and the whole rebuild costs minutes, so the
    standing rule is to regenerate rather than reason about whether a bump mattered. ``--check``
    is the cheap half: it reads version stamps and builds nothing, which is what CI runs.
    """
    from . import __version__, _build

    want = list(_build.TARGETS) if a.target == "all" else [a.target]

    if a.check:
        bad = _build.check(want)
        for tgt, f, got, exp in bad:
            print(f"{tgt}\t{f}\t{got}\t{exp}")
        say(f"{len(bad)} stale of {sum(len(_build.TARGETS[t][2]) for t in want)} artifact(s) "
            f"checked against {__version__}")
        raise SystemExit(1 if bad else 0)

    for tgt in want:
        label, fn, _files = _build.TARGETS[tgt]
        if fn is None:
            how = _build.EXTERNAL.get(tgt)
            say(f"{tgt}: {label} -- not buildable in this process; run: {how}" if how else
                f"{tgt}: {label} -- static, no generator on record; see data/PROVENANCE.md")
            continue
        with step(f"build {tgt}: {label}"):
            fn(say=lambda m: say(m, level=1, flush=True))
    say(f"shipped artifacts now stamped {__version__}; commit what changed")


def _cassette_rows(path, score_col, need_score=True):
    """Rows of a cassette / pool table, with the score column resolved and parsed.

    ``rank`` writes its aggregate under ``aggregate``; a hand-made table usually calls it ``score``.
    Rather than make the caller remember which, the default column is tried and the two known
    aliases are accepted, with the resolved name reported --- silently scoring the wrong column is
    the failure this avoids.
    """
    rows = _read_table(path)          # resolves `peptide` / `epitope`, and always sets `peptide`
    if not rows:
        raise SystemExit(f"{path}: no rows")
    cols = list(rows[0])
    # The restricting allotype, resolved and reported like the score column below. It is not
    # cosmetic: `select` keys its allotype coupling channel and its coverage on this column, so a
    # table spelling it `best_allele` -- every pipeline table does -- silently lands every unit on
    # one empty allotype, and the objective then prices no spread at all. `mm_allele_scored` leads
    # because that is the allele a `rank --prefix` row's numbers are actually against.
    if "allele" not in cols:
        acol = next((c for c in ("mm_allele_scored", "mm_allele", "best_allele", "allele_scored")
                     if c in cols), None)
        if acol:
            say(f"allele column: {acol!r}", level=1)
            for r in rows:
                r["allele"] = r.get(acol, "")
        else:
            say("no allele column: the allotype channel and coverage are off "
                f"(looked for allele, mm_allele_scored, mm_allele, best_allele; found {cols})",
                level=1)
    # The source gene, resolved the same way and for a reason that is not cosmetic either. It is
    # what `vector.self_origin_risk` excludes when it asks "does this register coincide with a
    # DIFFERENT expressed gene" -- so an empty one makes every register match its own source and the
    # safety screen withdraws the unit. Measured on a real donor whose table spelled it `gene_name`:
    # 18 of 20 units withdrawn over 20,150 findings, none of them a real off-target.
    if "gene" not in cols or not any(str(r.get("gene") or "").strip() for r in rows):
        gcol = next((c for c in ("mm_gene", "gene_name", "gene_symbol") if c in cols), None)
        if gcol:
            say(f"gene column: {gcol!r}", level=1)
            for r in rows:
                r["gene"] = r.get(gcol, "")
    # A pipeline table spells the peptide `epitope`, and one `rank --passthrough --prefix` has
    # annotated carries the caller's `epitope` beside our `mm_peptide` -- the same string. Resolved
    # like the score column below rather than made the caller's problem with a rename stage.
    if "peptide" not in cols:
        for r in rows:
            r["peptide"] = r.get("epitope", "")
    col = next((c for c in (score_col, "aggregate", "score", "epic") if c in cols), None)
    if col is None and need_score:
        raise SystemExit(f"{path}: no score column (looked for {score_col!r}, `aggregate`, `score`, "
                         f"`epic`; found {cols})")
    if col != score_col:
        say(f"score column: {col!r} (--score-column {score_col!r} not present)", level=1)
    for r in rows:
        try:
            r["_score"] = float(r[col])
        except (KeyError, TypeError, ValueError):
            r["_score"] = float("nan")
    return rows, col


def _f(v):
    """A table cell as a float, ``nan`` for an empty or unparseable one. ``nan`` is the right
    answer for a missing expression term -- :func:`mhcmatch.cassette.selectivity_delta` turns it
    into no preference rather than into a deleted candidate."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _by_donor(rows):
    """``{donor: [row]}`` preserving file order; one group called ``-`` when there is no column."""
    out = {}
    for r in rows:
        out.setdefault((r.get("donor") or r.get("patient") or "-").strip() or "-", []).append(r)
    return out


def cmd_cassette_select(a):
    """Choose k units from a donor's pool by maximising the mean-variance objective."""
    import numpy as np

    from . import cassette as CA

    rows, col = _cassette_rows(a.candidates, a.score_column)
    groups = _by_donor(rows)
    kw = {k: v for k, v in (("prevalence", a.prevalence), ("rho", a.rho), ("gamma", a.gamma))
          if v is not None}
    universe = [x.strip() for x in a.universe.split(",") if x.strip()] if a.universe else None
    say(f"{len(rows)} candidate(s) over {len(groups)} donor(s), scored on {col!r}; "
        f"k = {a.size} +/- {a.tol}", level=1)
    if a.block_live < 1.0:
        say(f"HLA loss priced at q = {a.block_live}: two units on one allotype are lost together",
            level=1)

    out, chosen = [], 0
    for donor, g in groups.items():
        bad = [i for i, r in enumerate(g) if not np.isfinite(r["_score"])]
        g = [r for i, r in enumerate(g) if i not in set(bad)]
        if bad:
            say(f"{donor}: {len(bad)} candidate(s) with no finite score, dropped", level=2)
        if not g:
            continue
        alleles = None if a.no_allele or "allele" not in g[0] else [r["allele"] for r in g]
        # `--floor` seeds every allotype the donor's own pool can supply; `--universe` is the true
        # genotype and is what lets `coverage` see an allotype holding zero units.
        uni = universe or (sorted(set(alleles)) if a.floor and alleles else None)
        sel_kw = dict(kw, block_live=a.block_live, universe=uni, max_share=a.max_share)
        if a.selectivity:
            missing = [c for c in ("expr_lvl", "expr_norm") if c not in g[0]]
            if missing:
                raise SystemExit(f"--selectivity needs the {' and '.join(missing)} column(s); "
                                 f"`mhcmatch rank` emits both. Found: {sorted(g[0])}")
            sel_kw.update(selectivity=a.selectivity,
                          expr_lvl=[_f(r["expr_lvl"]) for r in g],
                          expr_norm=[_f(r["expr_norm"]) for r in g])
        if a.no_dominance:
            sel_kw["dominance"] = False
        if a.rule == "v2":
            sel_kw.update(rule="v2", pi=a.not_worse, how=a.diversity)
        if a.feature_column:
            missing = [c for c in a.feature_column if c not in g[0]]
            if missing:
                raise SystemExit(
                    f"--feature-column names {' and '.join(missing)}, which the candidate table "
                    f"does not carry. `mhcmatch rank --score aggregate` emits expr_lvl, expr_norm, "
                    f"C_phys_buried and C_phys_charge. Found: {sorted(g[0])}")
            sel_kw.update(features=np.array([[_f(r[c]) for c in a.feature_column] for r in g]),
                          feature_names=tuple(a.feature_column))
        if a.coexpr_gtex:
            if "gene" not in g[0]:
                raise SystemExit("--coexpr-gtex needs the `gene` column; `mhcmatch rank` emits it")
            from . import expression as EX
            genes = [r.get("gene", "") for r in g]
            sel_kw["coexpr"] = EX.coexpression(genes)
            known = sum(1 for x in genes if x)
            say(f"{donor}: GTEx co-expression over {known} of {len(genes)} unit(s) with a gene "
                "symbol; the rest contribute no pair information", level=2)
        size, tol = a.size, a.tol
        if a.confidence is not None:
            # `-k` becomes the manufacturing ceiling and the donor's own pool sets the size. A
            # donor whose head of list is weak needs more units to reach the same confidence, and
            # the ceiling is where that answer gets reported rather than silently rounded down.
            # `block_live` too: the probe asks how many units reach a confidence, and a cassette
            # that can lose a whole allotype at once needs more of them. Without it `--confidence`
            # sized every donor as if `--block-live` were 1.0 while `select` below honoured it, so
            # the two flags did not compose and the size came out low.
            probe = CA.size_for([r["_score"] for r in g], [r["peptide"] for r in g], alleles,
                                target=a.target, confidence=a.confidence, k_max=a.size,
                                block_live=(a.block_live if alleles is not None else 1.0),
                                **{k: v for k, v in kw.items() if k != "gamma"})
            size, tol = probe["k"], 0
            say(f"{donor}: {size} unit(s) for P(>= {a.target}) >= {a.confidence:.2f}"
                + ("" if probe["reached"] else
                   f" -- NOT REACHED, capped at k = {a.size}, best {probe['p_at_least']:.3f}"),
                level=1)
        c = CA.select([r["_score"] for r in g], [r["peptide"] for r in g], alleles,
                      k=size, tol=tol, **sel_kw)
        if a.selectivity:
            # The counterfactual is what makes a stated weight auditable rather than a knob: what
            # the same pool would have built at w = 0, and what the trade cost in expected units.
            base = CA.select([r["_score"] for r in g], [r["peptide"] for r in g], alleles,
                             k=size, tol=tol, **dict(sel_kw, selectivity=0.0))
            d = CA.selectivity_delta([_f(r["expr_lvl"]) for r in g],
                                     [_f(r["expr_norm"]) for r in g])
            say(f"{donor}: selectivity w={a.selectivity} traded yield {base.yield_:.3f} -> "
                f"{c.yield_:.3f} unit(s) for mean tumour-over-normal "
                f"{d[base.index].mean():+.3f} -> {d[c.index].mean():+.3f} log2-fold; "
                f"{len(set(c.index) - set(base.index))} of {c.k} slot(s) changed", level=1)
        if c.trimmed:
            say(f"{donor}: pool of {c.pool_n} trimmed to {c.pool_n - c.trimmed} before the "
                f"coupling matrix (see mhcmatch.cassette.MAX_POOL)", level=1)
        chosen += c.k
        for slot, (i, pi) in enumerate(zip(c.index, c.p), start=1):
            r = dict(g[i])
            r.pop("_score", None)
            # --passthrough: the caller's own columns lead, so the chosen units carry whatever the
            # candidate table carried -- the long window a cassette is actually built from
            # (`epitope_context`), the variant class the quota reads, the caller's identifiers.
            # Without it those have to be joined back on the peptide, and a pool may hold the same
            # peptide on two allotypes.
            carried = r if getattr(a, "passthrough", False) else {}
            out.append({**carried,
                        "donor": donor, "slot": slot, "peptide": r["peptide"],
                        "allele": r.get("allele", ""), "gene": r.get("gene", ""),
                        "score": f"{g[i]['_score']:.6f}", "p": f"{pi:.6f}",
                        "k": c.k, "pool_n": c.pool_n, "offset": f"{c.offset:.6f}",
                        "energy": f"{c.energy:.6f}", "lam": f"{c.lam:.6f}",
                        "rho": c.rho, "gamma": c.gamma, "channels": "+".join(c.channels),
                        "block_live": c.block_live, "selectivity": c.selectivity,
                        "rule": c.rule, "pi": c.pi, "not_worse": f"{c.not_worse:.4f}",
                        "diversity": f"{c.diversity:.4f}",
                        "n_covered": c.coverage.get("n_covered", ""),
                        "n_allotypes": c.coverage.get("n_allotypes", "")})
        # A tolerance that is spent is a result, not a detail: the objective has an internal
        # optimum size, and it moves with the prevalence and with rho. Saying so is cheaper than
        # letting somebody discover that `-k 20 --tol 5` returned 15 and wonder whether it broke.
        if c.k != size:
            say(f"{donor}: {c.k} units, not {size} -- the objective peaks there inside the "
                f"tolerance (adding the next unit costs more variance than it buys in mean)",
                level=1)
        if c.rule == "v2":
            say(f"{donor}: v2 moved {c.swaps} of {c.k} slot(s) off the sort for diversity "
                f"{c.diversity:.4f} ({c.how}), keeping P(not worse) = {c.not_worse:.3f} "
                f"against the stated floor {c.pi:.2f}", level=1)
        cov = c.coverage
        say(f"{donor}: {c.k} of {c.pool_n}, yield {c.yield_:.3f} unit(s), lam {c.lam:+.3f} nats"
            + (f", {cov['n_covered']}/{cov['n_allotypes']} allotype(s) covered, "
               f"{cov['entropy_ratio']:.3f} of maximum entropy" if cov else ""), level=2)
    say(f"selected {chosen} unit(s) over {len(groups)} donor(s)", level=1)
    _write_rows(out, a.out)


def cmd_cassette_score(a):
    """Score finished cassettes on axes that survive changing donor and changing size."""
    from . import cassette as CA

    rows, col = _cassette_rows(a.cassettes, a.score_column)
    cass = _by_donor(rows)
    pool = _by_donor(_cassette_rows(a.pool, a.score_column)[0]) if a.pool else {}
    prev = POOL_PREVALENCE if a.prevalence is None else a.prevalence

    # The offset is the whole question. One over the file makes two cassettes' `yield` comparable
    # levels; one per donor makes it an enrichment against that donor's own background, and no two
    # donors' numbers are then on one axis. Neither is wrong; the flag is which one was meant.
    if a.per_donor_offset:
        offs = {d: CA.prob_offset([r["_score"] for r in (pool.get(d) or g)], prev)
                for d, g in cass.items()}
        say(f"one offset per donor over {len(offs)} donor(s): every donor's mean probability is "
            f"now {prev:.4f} by construction, so `yield` is an ENRICHMENT, not a level", level=1)
    else:
        src = [r["_score"] for g in (pool or cass).values() for r in g]
        b = CA.prob_offset(src, prev)
        offs = {d: b for d in cass}
        say(f"one offset over {len(src)} row(s) of the {'pool' if pool else 'cassette'} file: "
            f"b = {b:+.4f} at prevalence {prev:.4f}", level=1)

    out = []
    for donor, g in cass.items():
        pg = pool.get(donor)
        if a.pool and pg is None:
            raise SystemExit(f"--pool has no rows for donor {donor!r}; every scored cassette needs "
                             "the pool it was chosen from, or drop --pool")
        alleles = [r["allele"] for r in g] if "allele" in g[0] else None
        kw = {"rho": a.rho} if a.rho is not None else {}
        s = CA.score([r["_score"] for r in g], [r["peptide"] for r in g], alleles,
                     pool_scores=None if pg is None else [r["_score"] for r in pg],
                     pool_peptides=None if pg is None else [r["peptide"] for r in pg],
                     offset=offs[donor], block_live=a.block_live, target=a.target,
                     universe=[x.strip() for x in a.universe.split(",") if x.strip()]
                     if a.universe else None, **kw)
        cov = s.pop("coverage", {}) or {}
        s = {"donor": donor, **s,
             "n_allotypes": cov.get("n_allotypes", ""), "allotype_gini": cov.get("gini", ""),
             "allotype_entropy_ratio": cov.get("entropy_ratio", "")}
        out.append({k: (f"{v:.6f}" if isinstance(v, float) else ("" if v is None else v))
                    for k, v in s.items()})
    say(f"scored {len(out)} cassette(s), sizes {min(len(g) for g in cass.values())}"
        f"-{max(len(g) for g in cass.values())}", level=1)
    _write_rows(out, a.out)


def _write_rows(rows, out):
    """One TSV, header from the first row, to ``out`` or stdout. Empty input writes nothing."""
    if not rows:
        say("nothing to write", level=1)
        return
    cols = list(rows[0])
    fh = open(out, "w") if out else sys.stdout
    try:
        print("\t".join(cols), file=fh)
        for r in rows:
            print("\t".join(str(r.get(c, "")) for c in cols), file=fh)
    finally:
        if out:
            fh.close()
            say(f"wrote {len(rows)} row(s) -> {out}", level=1)


def _add_deslip_opts(p) -> None:
    """The `deslip` options, on one parser -- shared by `cassette deslip` and the alias."""
    p.add_argument("cds", help="coding sequence, or a FASTA path (T or U, case-insensitive)")
    p.add_argument("--fix", metavar="FILE",
                   help="write the repaired CDS here: every TTT before a T/C-starting codon becomes "
                        "TTC, which is synonymous, so the protein is unchanged")
    p.add_argument("--out", metavar="FILE", help="write the site TSV here instead of stdout")


def _add_vector_opts(p, require_n0: bool = True) -> None:
    """Every option of the assembly pipeline, on one parser.

    Called for ``cassette build``, for ``cassette order`` (which does not select, so
    ``--n0`` is not required), and for the deprecated ``vector`` alias. One definition, so
    the three cannot drift into three different flag sets.
    """
    p.add_argument("--candidates", required=True, metavar="FILE",
                    help="TSV of units: peptide, gene, allele, p (+ optional mutation_index, cls). "
                         "`peptide` is the LONG window around the mutation, not the minimal epitope "
                         "-- a minimal peptide loads onto any cell without costimulation and is the "
                         "tolerising configuration. `-` = stdin")
    p.add_argument("--context", metavar="FILE",
                    help="the window FASTA `rank` was run on. With it, --candidates may be `rank`'s "
                         "own output of MINIMAL epitopes: each is joined back to its source window "
                         "and re-centred as a long unit, one per variant rather than one per "
                         "register. Without it --candidates must already carry long windows")
    p.add_argument("--unit-column", metavar="COL",
                    help="the column of --candidates holding the LONG window, when the table "
                         "carries it beside the minimal epitope (a pipeline table spells it "
                         "`epitope_context`, and at 27 residues it is already --unit-length). The "
                         "alternative to --context, which rebuilds the window from the variant "
                         "FASTA when the table has none; ignored when --context is given")
    p.add_argument("--unit-length", type=int, default=27, metavar="N",
                    help="unit window length for --context (default 27, the BioNTech backbone "
                         "configuration; see mhcmatch.vector.unit)")
    p.add_argument("--n0", type=float, required=require_n0, metavar="F",
                    help="per-allotype capacity, the one free parameter of the stopping rule. "
                         "REQUIRED and with no default on purpose: nothing in the public record fits "
                         "it, so the value is yours to defend and it is recorded in the output")
    p.add_argument("--quota", metavar="ARM=SLOTS:TARGET",
                    help="compose the cassette to quotas instead of taking the ranked top, e.g. "
                         "'mhc1=8:2,mhc2=4:1,nonconventional=3:1' -- eight class-I slots of which "
                         "at least two should respond, and so on. Arms are disjoint: a unit whose "
                         "`kind` column is anything but `missense` is charged to `nonconventional`, "
                         "so the constraint bites. The same slot budgets filled by score alone are "
                         "reported beside it, because 'not the same as ranking' is a claim about "
                         "YOUR candidates")
    p.add_argument("--block-live", type=float, default=0.5, metavar="Q",
                    help="P(a block is live) in the response model behind --quota (default "
                         "%(default)s). A block is an allotype: if the recipient never mounts a "
                         "response on that allotype, none of its units respond however good they "
                         "are. Measure it on your own readout with "
                         "mhcmatch.portfolio.betabinom_rho before trusting a default")
    p.add_argument("--evenness", type=float, default=0.0, metavar="W",
                    help="weight on class-I allotype evenness (H/Hmax) in --quota's objective "
                         "(default %(default)s = off). The block model already prefers spread when "
                         "spread helps; this is for when it does not and you want it anyway. "
                         "Homozygosity is handled: the denominator is the DISTINCT allotypes in "
                         "--alleles, so a homozygous locus is not scored as a design flaw")
    p.add_argument("--alleles", help="the recipient's allotypes for junction scoring "
                                      "(comma-separated or a file); default = those in the table")
    p.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    p.add_argument("--cls-filter", action="store_true",
                    help="select only units whose own `cls` matches --cls")
    p.add_argument("--screen", action="store_true",
                    help="withdraw units on essential-tissue risk BEFORE selecting. Costs a "
                         "whole-proteome index (minutes, several GB); without it no safety check "
                         "runs at all and the cassette carries whatever it was handed")
    p.add_argument("--screen-mode", default="veto", choices=("veto", "graded"),
                    help="what --screen does with a finding (default %(default)s). `veto` withdraws "
                         "the unit, as shipped. `graded` withdraws only findings at or above "
                         "--veto-tpm and keeps the rest as a per-unit OFF-TARGET FINGERPRINT, "
                         "reported in the reason rows and priced into --quota by "
                         "--weight-offtarget. A 27-mer carries ~70 registers, so letting any one of "
                         "them veto the unit is not the specificity the per-register measurement "
                         "reads as")
    p.add_argument("--min-tpm", type=float, default=0.25, metavar="F",
                    help="essential-tissue expression floor for --screen: the REPORTING floor, "
                         "below which a finding is not recorded at all. 0.25 because MAGE-A12 sits "
                         "at 0.33 TPM in brain and killed two patients; a conventional 5 would pass "
                         "it")
    p.add_argument("--veto-tpm", type=float, default=5.0, metavar="F",
                    help="essential-tissue expression at or above which a finding WITHDRAWS the "
                         "unit under --screen-mode graded (default %(default)s, the conventional "
                         "'is it expressed' cut). Distinct from --min-tpm on purpose: 0.25 TPM is "
                         "'detectable somewhere', which nearly every human gene is, so it is a "
                         "floor for reporting and not a line for exclusion")
    p.add_argument("--weight-offtarget", type=float, default=0.0, metavar="W",
                    help="price of one off-target fingerprint entry in --quota's objective "
                         "(default %(default)s = off). The composed value becomes "
                         "P(X >= target) - W * sum(distinct essential-tissue genes a unit reaches). "
                         "Charged to the objective, never to the unit's calibrated p")
    p.add_argument("--max-subs", type=int, default=0, metavar="N",
                    help="self-origin search radius for --screen. 0 = exact coincidence, which is "
                         "the default because the decision is per unit while the search is per "
                         "register: at radius 1 over 8-11mers, 3 of 6 random 27-mers get withdrawn "
                         "by chance. Raise it only together with dropping 8-mers")
    p.add_argument("--report-subs", type=int, default=0, metavar="D", choices=(0, 1),
                    help="also REPORT (never withdraw) registers within D substitutions of a "
                         "different, non-homologous, expressed gene (default %(default)s = off; "
                         "1 is the only other value). The two TCR deaths were near-identity, not "
                         "identity, but a radius-1 veto costs two thirds of every cassette to buy "
                         "a hazard --max-subs 0 has largely already taken: measured on 178 "
                         "validated immunogenic neoantigens, exact withdraws 1.1%% and d=1 to any "
                         "expressed gene reaches 70.2%%. So d=1 findings ride along in the "
                         "fingerprint instead, priced only by --weight-offtarget")
    p.add_argument("--report-identity", type=float, default=0.5, metavar="F",
                    help="flanking-identity ceiling for --report-subs (default %(default)s). A d=1 "
                         "match to a gene that shares the unit's flanks is descent, not mimicry, "
                         "and tolerance already covers it; this cut is what takes the reported "
                         "tier's different-gene hits from 230 to 74 at L=9, 130 of the 156 it "
                         "removes being one locus under two symbols. It separates loci, not "
                         "superfamilies: a 27-mer bounds the comparison at +/-9-10 residues, so "
                         "NRAS -> KRAS survives at 0.23 and is reported")
    p.add_argument("--report-threshold", type=float, default=-1.47712, metavar="F",
                    help="-log10(%%rank) the off-target variant must itself reach to be reported "
                         "under --report-subs (default %(default)s = 30%% rank). A d=1 coincidence "
                         "no allotype presents is a sequence coincidence, not a safety "
                         "consideration. The cut looks permissive because it is read off the "
                         "positives: 97.2%% of 176 assayed immunogenic neoantigens clear 30%% on "
                         "this scorer, where the conventional 2%% would discard three in ten of "
                         "them -- the wrong error to make on a safety question")
    p.add_argument("--objective", default="sum", choices=("sum", "rate"),
                    help="junction cost: `sum` of the strongest binder per junction (pVACvector's "
                         "logic, biased toward the shortest spacer), or `rate` = binders per "
                         "register, which is length-neutral and needs --binder-threshold. The two "
                         "disagree on real payloads, so choose")
    p.add_argument("--binder-threshold", type=float, metavar="F",
                    help="-log10(%%rank) above which a junction window counts as a binder; "
                         "required by --objective rate")
    p.add_argument("--threshold", type=float, metavar="F",
                    help="stop at the first spacer whose worst junction falls at or below this, "
                         "instead of trying them all and taking the cheapest")
    p.add_argument("--linker", metavar="NAME",
                    help="PIN one linker instead of sweeping: a preset name (`mhcmatch cassette "
                         "linkers` lists them), explicit residues, or `none` for direct "
                         "concatenation. Use this when the construct format is already decided "
                         "and only the ordering is open; without it every spacer is tried and the "
                         "cheapest wins")
    p.add_argument("--mrna", metavar="FILE",
                    help="also write the assembled mRNA as FASTA: 5' UTR, start codon, leader, the "
                         "cassette CDS, trailer, stop, 3' UTR and poly(A), in that order, with a "
                         "parts map on the header line. Every backbone element defaults to nothing "
                         "and is supplied by the flags below -- this library does not invent a UTR")
    p.add_argument("--leader", metavar="SEQ", default="",
                    help="amino acids translated in frame BEFORE the cassette (a secretory signal "
                         "peptide); a sequence or a FASTA/plain file")
    p.add_argument("--trailer", metavar="SEQ", default="",
                    help="amino acids translated in frame AFTER the cassette (an MHC class-I "
                         "trafficking domain); a sequence or a FASTA/plain file")
    p.add_argument("--utr5", metavar="SEQ", default="",
                    help="5' untranslated region, NUCLEOTIDES; a sequence or a file")
    p.add_argument("--utr3", metavar="SEQ", default="",
                    help="3' untranslated region, NUCLEOTIDES; a sequence or a file")
    p.add_argument("--poly-a", type=int, default=0, metavar="N",
                    help="length of the template-encoded poly(A) tail in adenosines (default 0, "
                         "i.e. none -- the tail is a property of the vector, not of the payload)")
    p.add_argument("--fasta", metavar="FILE", help="also write the cassette sequence as FASTA")
    p.add_argument("--fasta-nt", metavar="FILE",
                    help="also write the cassette CODING SEQUENCE as FASTA -- highest-usage human "
                         "codon per residue, backed off to avoid homopolymers, then deslipped. "
                         "Epitope cassette only: no start, no stop, no leader, no trafficking "
                         "domain")
    p.add_argument("--map", metavar="FILE", dest="map_tsv",
                    help="also write the cassette MAP as TSV: one row per unit, linker and "
                         "predicted epitope, with 1-based coordinates over the cassette, the "
                         "presenting allele, and which class-I and class-II epitopes overlap each "
                         "other. A peptide presented by two of the recipient's alleles gets TWO "
                         "rows -- at a heterozygous locus those are two presentation events")
    p.add_argument("--map-json", metavar="FILE",
                    help="the same map as JSON, plus the per-unit summary and the sequence, which "
                         "is what a viewer needs to draw the cassette without recomputing anything")
    p.add_argument("--map-binder", choices=("strong", "weak"), default="weak",
                    help="which NetMHCpan cut-off the map annotates. NetMHCpan calls class I strong "
                         "at %%rank <= 0.5 and weak at <= 2.0; NetMHCIIpan calls class II strong at "
                         "<= 2.0 and weak at <= 10.0 -- so the two classes do NOT share a number. "
                         "Default `weak`, because the map reports and never selects (default: weak)")
    p.add_argument("--map-threshold", type=float, default=None, metavar="F",
                    help="override the class-I %%rank cut-off for the map (default: from "
                         "--map-binder, so 2.0 weak / 0.5 strong)")
    p.add_argument("--map-threshold-mhc2", type=float, default=None, metavar="F",
                    help="override the class-II %%rank cut-off (default: from --map-binder, so "
                         "10.0 weak / 2.0 strong)")
    p.add_argument("--map-alleles-mhc2", metavar="LIST",
                    help="the recipient's class-II allotypes (comma-separated or a file). Without "
                         "them the map carries class I only, and a unit's `self_help` column -- "
                         "whether its CD8 epitope has overlapping CD4 help from the SAME unit -- "
                         "cannot be computed")
    _add_store_opts(p)
    p.add_argument("--out", metavar="FILE", help="write the report TSV here instead of stdout")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mhcmatch", description="peptide-MHC presentation tools")
    # `bench/run_epic.sh` gates the whole reproduction on this matching the checkout's pyproject,
    # so it has to exist and it has to print the installed distribution's version, not a constant.
    ap.add_argument("--version", action="version", version=f"mhcmatch {__version__}")
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
    s.add_argument("--out", help="write TSV here instead of the aligned text")
    s.add_argument("--tsv", action="store_true",
                    help="emit TSV to stdout rather than the aligned text")
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
    lg.add_argument("--out", help="write TSV here instead of the aligned text")
    lg.add_argument("--tsv", action="store_true",
                    help="emit TSV to stdout rather than the aligned text")
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
    pr.add_argument("--background", default="proteome", choices=("ligand", "ligand-pooled", "proteome", "markov"))
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

    al = sub.add_parser("alleles",
                        help="a donor's HLA typing file -> the allele list --alleles accepts")
    al.add_argument("input", help="typing TSV with an `Allele` column (OptiType / kourami / HLA-LA "
                                  "and a donor's own .alleles.tsv), a comma-separated list, or one "
                                  "name per line")
    al.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"),
                    help="which class to emit. One typing file holds both, and a class-I panel "
                         "handed a DQB1 name resolves it to nothing (default: %(default)s)")
    al.add_argument("--form", default="key", choices=("key", "input"),
                    help="`key` (default) emits the pseudosequence key every scoring path uses; "
                         "`input` emits the typed name that resolved to it, for a caller that wants "
                         "its own spelling back")
    al.add_argument("--out", help="write the comma-separated list here instead of stdout")
    al.set_defaults(fn=cmd_alleles)

    rk = sub.add_parser("rank", help="rank neoantigen candidates (FASTA of windows, or a scored table)")
    rk.add_argument("mode", nargs="?", choices=("fasta", "table", "pairs"),
                    help="fasta: mutation-spanning window FASTA + donor alleles. "
                         "table: a .scored.csv already produced by another tool. "
                         "pairs: a TSV of peptide / wt_peptide / allele (+ optional gene, tpm) -- "
                         "the shape a neoantigen screen is distributed in. "
                         "Omit with --coefficients / --holdout")
    rk.add_argument("input", nargs="?",
                    help="the .peptide.fasta, the .scored.csv, or the pairs TSV (.gz ok, - = stdin)")
    rk.add_argument("--coefficients", action="store_true",
                    help="print the fitted aggregate as TSV (block, term, coef, sd, boot_sd, z, "
                         "p, 95%% interval, sign stability) and score nothing. Read from the "
                         "shipped data/aggregate_mhc1.json, never refitted")
    rk.add_argument("--holdout", action="store_true",
                    help="print that fit's leave-one-screen-out and cross-validated AUROCs as "
                         "TSV and score nothing")
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
    rk.add_argument("--passthrough", action="store_true",
                    help="mode=pairs/table: emit every column of the input table, unchanged and in "
                         "its own order, ahead of this command's own -- so a caller's table comes "
                         "back annotated and re-ordered by the aggregate rather than replaced. Do "
                         "NOT do this with a join instead: a cell naming several alleles is split "
                         "and the best presenter stands for the row, so the output shares neither "
                         "its length nor its allele column with the input. Rows come out in `rank` "
                         "order, not the listing order an exact known-epitope match floats to")
    rk.add_argument("--prefix", default="", metavar="STR",
                    help="with --passthrough, prefix the columns THIS command adds (`mm_` is what "
                         "the shipped deliverables use), so a name it shares with the caller's "
                         "table -- `score`, `allele`, `rank` -- does not appear twice")
    rk.add_argument("--context", metavar="FILE",
                    help="mode=pairs: the window FASTA the candidates were called on, read for the "
                         "germline arm of each window (`wt_window`) so agretopicity and "
                         "`d_occupancy` are defined. A candidate table usually carries the mutant "
                         "k-mer and nothing the wild type is recoverable from; without this every "
                         "row is `wt_absent`. Only equal-length window pairs are used, so a "
                         "frameshift, a fusion and an indel stay wild-type-less, which is what "
                         "they are")
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

    # ---------------------------------------------------------------- cassette
    # The one command group with sub-verbs. Everything about a cassette lives under it: what goes
    # in (`select`), what a finished one is worth (`score`), and how it is assembled (`build`,
    # `order`, `deslip`). `vector` and `deslip` remain as deprecated top-level aliases for one
    # release, because they are in every published pipeline config we know of.
    ca = sub.add_parser("cassette",
                        help="choose the units of a vaccine cassette, score one, and assemble it")
    casub = ca.add_subparsers(dest="cassette_cmd", required=True)

    cs = casub.add_parser("select",
                          help="choose k units (+/- tol) from a donor's candidate pool, maximising "
                               "the mean-variance objective rather than sorting on the score")
    cs.add_argument("--candidates", required=True, metavar="FILE",
                    help="TSV of the donor's WHOLE candidate pool: peptide, score (+ optional "
                         "allele, gene, id). `mhcmatch rank`'s output is the intended input. `-` = "
                         "stdin. Do not pass a shortlist already filtered on binding or expression: "
                         "those carry the two largest coefficients in the model, so a pool cut on "
                         "them has no range left along them")
    cs.add_argument("-k", "--size", type=int, default=20, metavar="N",
                    help="target cassette size (default: 20)")
    cs.add_argument("--tol", type=int, default=0, metavar="N",
                    help="manufacturing tolerance: the reported size is the one in [k-tol, k+tol] "
                         "with the largest objective, with the lower bound raised to the number of "
                         "--universe allotypes the pool can supply (default: 0, i.e. exactly k)"
                         "with the largest objective (default: 0, i.e. exactly k)")
    cs.add_argument("--prevalence", type=float, metavar="P",
                    help=f"pool response prevalence the probabilities are anchored on "
                         f"(default: {POOL_PREVALENCE:.4f}, TESLA's 37 of 615). Fitted ONCE over "
                         "the pool and held; halving it roughly halves every probability and moves "
                         "no rank")
    cs.add_argument("--rho", type=float, metavar="R",
                    help=f"intra-cassette response correlation (default: {CA_RHO}, IVAC MUTANOME). "
                         "Measure your own with mhcmatch.portfolio.betabinom_rho")
    cs.add_argument("--gamma", type=float, metavar="G",
                    help="risk aversion: one unit of variance in the responding-unit count is worth "
                         "this many expected units, PER UNIT of the cassette. A stated design "
                         "preference; the default 1.0 is divided by the design effect 1+rho(k-1) "
                         "so it means the same trade at every k, and passing this uses the number "
                         "given instead")
    cs.add_argument("--confidence", type=float, metavar="C",
                    help="size each cassette to the donor instead of using -k: the smallest number "
                         "of units reaching P(>= --target responses) >= C for that donor's own "
                         "pool, with -k as the manufacturing ceiling. A donor whose best "
                         "candidates are weak needs more units to reach the same C, and is "
                         "reported at the ceiling when C is out of reach")
    cs.add_argument("--target", type=int, default=1, metavar="M",
                    help="how many responding units --confidence is about (default: 1)")
    cs.add_argument("--block-live", type=float, default=1.0, metavar="Q",
                    help="how often each allotype survives -- the HLA-loss rate. Below 1 the "
                         "objective prices losing an allele: two units on one allotype are lost "
                         "together, so the coupling between them is the exact covariance that "
                         "implies. A unit whose p exceeds q raises rather than being clipped "
                         "(default: 1.0, nothing is ever lost)")
    cs.add_argument("--floor", action="store_true",
                    help="give every allotype the donor's pool can supply one unit before the free "
                         "slots are filled. A manufacturing constraint, not an objective term")
    cs.add_argument("--max-share", type=float, metavar="F",
                    help="no allotype may hold more than this share of the cassette (0.4 at k=20 "
                         "caps each at 8 units). Refuses rather than relaxing when the share and "
                         "the floor cannot both hold")
    cs.add_argument("--universe", metavar="LIST",
                    help="the donor's DISTINCT allotypes, comma-separated -- the denominator "
                         "coverage is reported against, so an allotype holding zero units is "
                         "visible. Also sets the floor. Without it, coverage is taken over the "
                         "labels the cassette happens to carry and cannot see one it missed")
    cs.add_argument("--selectivity", type=float, default=0.0, metavar="W",
                    help="expected responding units per log2-fold of tumour-over-normal abundance, "
                         "charged to the objective as w*(expr_lvl - expr_norm) and NEVER to p. A "
                         "stated design preference like --gamma: the shipped model fits both "
                         "expression terms POSITIVE because it was fitted on `will this respond`, "
                         "and `high in tumour, low in normal` is a different question (default: 0)")
    cs.add_argument("--rule", choices=("v1", "v2"), default="v1",
                    help="v1 maximises the mean-variance objective. v2 takes the plain sort as a "
                         "floor and spends the slack on diversity: because p is a probability, "
                         "many size-k sets are indistinguishable in how many units respond, and v2 "
                         "returns the most diverse set that is still, with probability at least "
                         "--not-worse, no worse than sorting (default: v1)")
    cs.add_argument("--not-worse", type=float, default=0.5, metavar="P",
                    help="v2 only: the floor on P(this cassette catches at least as much as the "
                         "top-k sort). 1.0 returns the sort itself; lower buys more diversity and "
                         "says how often you are willing to be wrong (default: 0.5)")
    cs.add_argument("--diversity", choices=("minmax", "mean"), default="minmax",
                    help="v2 only: how diversity is aggregated over the four axes. `minmax` "
                         "maximises the WORST-covered axis, since a cassette is undone by its worst "
                         "shared failure mode; `mean` averages them, which dilutes any single axis "
                         "that separates (default: minmax)")
    cs.add_argument("--feature-column", action="append", metavar="COL",
                    help="add a coupling channel on this per-unit column, so two units alike on it "
                         "cost each other. Repeatable. `rank --score aggregate` emits the four "
                         "worth trying: C_phys_buried, C_phys_charge, expr_lvl, expr_norm. Unlike "
                         "--selectivity this changes what a PAIR costs, not what a unit is worth "
                         "(default: none, and the objective sees only the score)")
    cs.add_argument("--coexpr-gtex", action="store_true",
                    help="couple two units whose source genes share a GTEx tissue profile, so a "
                         "cassette does not spend two slots on one transcriptional programme. "
                         "Reads the `gene` column; a gene the panel does not carry contributes no "
                         "pair information rather than being dropped")
    cs.add_argument("--no-dominance", action="store_true",
                    help="drop the score-dominance channel, leaving only mechanism-based ones. It "
                         "is the one channel built from the score rather than from biology, and "
                         "the statistic it corresponds to fits ATTRACTIVE on the observational "
                         "arm, where the greedy 1-1/e guarantee does not hold")
    cs.add_argument("--score-column", default="score", metavar="COL",
                    help="which column holds the aggregate log-odds (default: score, which is "
                         "what `rank` writes; for a `rank --passthrough --prefix mm_` table "
                         "pass `mm_score`)")
    cs.add_argument("--no-allele", action="store_true",
                    help="ignore the allele column, so the overlap has no allotype channel. What a "
                         "trial that published no per-patient genotype is left with")
    cs.add_argument("--passthrough", action="store_true",
                    help="emit every column of --candidates ahead of this command's own, so the "
                         "chosen units keep the long window a cassette is built from and the "
                         "variant class `--quota` reads. This command's own columns win a name "
                         "clash, because they are the ones downstream reads")
    cs.add_argument("--out", metavar="FILE", help="write the chosen units here instead of stdout")
    cs.set_defaults(fn=cmd_cassette_select)

    cq = casub.add_parser("score",
                          help="score finished cassettes -- across donors and across sizes -- on "
                               "expected responding units, P(>=1) under the block model, and lambda")
    cq.add_argument("--cassettes", required=True, metavar="FILE",
                    help="TSV of manufactured units: peptide, score (+ optional donor, allele). "
                         "Rows are grouped by `donor` when the column is present, so one file may "
                         "hold many cassettes of different sizes. `-` = stdin")
    cq.add_argument("--pool", metavar="FILE",
                    help="the candidate pool each cassette was chosen from, same columns. With it "
                         "you also get `lam`: nats above a uniform random subset of that donor's own "
                         "pool, which is the axis that survives changing donor AND changing size")
    cq.add_argument("--prevalence", type=float, metavar="P",
                    help="prevalence the shared offset is anchored on (see `cassette select`)")
    cq.add_argument("--per-donor-offset", action="store_true",
                    help="fit one offset per donor instead of one over the whole file. This reports "
                         "an ENRICHMENT, not a level: every donor's mean probability becomes the "
                         "prevalence by construction, so the numbers are no longer probabilities "
                         "and no two donors' offsets are comparable")
    cq.add_argument("--rho", type=float, metavar="R", help="intra-cassette response correlation")
    cq.add_argument("--block-live", type=float, default=1.0, metavar="Q",
                    help="how often each block is live. A unit whose marginal p exceeds it raises "
                         "rather than being clipped (default: 1.0)")
    cq.add_argument("--universe", metavar="LIST",
                    help="the donor's DISTINCT allotypes, comma-separated. Without it `coverage` is "
                         "computed over the labels the cassette carries, so an allotype holding "
                         "zero units is invisible -- which is the inequality the index exists for")
    cq.add_argument("--target", type=int, default=1, metavar="M",
                    help="report P(at least M responses) (default: 1)")
    cq.add_argument("--score-column", default="score", metavar="COL")
    cq.add_argument("--out", metavar="FILE", help="write the report TSV here instead of stdout")
    cq.set_defaults(fn=cmd_cassette_score)

    cb = casub.add_parser("build",
                          help="assemble a cassette end to end: withdraw on safety, size each "
                               "allotype, order the units, pick a spacer, emit a map")
    _add_vector_opts(cb)
    cb.set_defaults(fn=cmd_vector)

    co = casub.add_parser("order",
                          help="the assembly half alone: order units already chosen and pick the "
                               "spacer, minimising junctional binding")
    _add_vector_opts(co, require_n0=False)
    co.set_defaults(fn=cmd_vector, order_only=True)

    cl = casub.add_parser("linkers",
                          help="list the named linker presets --linker accepts")
    cl.set_defaults(fn=cmd_linkers)

    cd = casub.add_parser("deslip",
                          help="find (and repair) the m1-pseudouridine +1 frameshift motif in a "
                               "cassette coding sequence")
    _add_deslip_opts(cd)
    cd.set_defaults(fn=cmd_deslip)

    # ---------------------------------------------------------------- deprecated aliases
    vc = sub.add_parser("vector", help="DEPRECATED alias for `cassette build`")
    _add_vector_opts(vc)
    vc.set_defaults(fn=cmd_vector, deprecated="cassette build")

    ds = sub.add_parser("deslip", help="DEPRECATED alias for `cassette deslip`")
    _add_deslip_opts(ds)
    ds.set_defaults(fn=cmd_deslip, deprecated="cassette deslip")

    gn = sub.add_parser("genes",
                        help="add a `gene` column to a peptide table -- the parent gene each "
                             "candidate derives from, which is what the expression terms are "
                             "keyed on")
    gn.add_argument("input", help="TSV(.gz) with a peptide column (`-` = stdin); every other "
                                  "column is carried through unchanged")
    gn.add_argument("--out", metavar="FILE", help="write TSV here instead of stdout")
    gn.add_argument("--peptide-col", default="peptide", help="the column holding the peptide")
    gn.add_argument("--species", default="human",
                    help="reference proteome: human / mouse / a pathogen stem auto-fetched from "
                         "the public dataset, or a FASTA(.gz) path")
    gn.add_argument("--max-subs", type=int, default=2,
                    help="search radius. 2 by default because a neoantigen can carry more than "
                         "one mutation; at radius 1 one screen resolved 88.2%% of its peptides "
                         "against 96.8%% at radius 2. Only the nearest shell votes")
    _add_thread_opt(gn)
    gn.set_defaults(fn=cmd_genes)

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
    xp.add_argument("--out", help="write TSV here instead of the aligned text")
    xp.add_argument("--tsv", action="store_true",
                    help="emit TSV to stdout rather than the aligned text")
    xp.set_defaults(fn=cmd_expression)

    bl = sub.add_parser("build",
                        help="regenerate the shipped data artifacts (release task), or --check them")
    bl.add_argument("target", nargs="?", default="all",
                    choices=("all", "anchor", "corpus", "recognition"),
                    help="what to rebuild (default: all)")
    bl.add_argument("--check", action="store_true",
                    help="do not build: report artifacts whose version stamp is behind "
                         "__version__, one per line, and exit 1 if any are. What CI runs")
    bl.set_defaults(fn=cmd_build)

    # Every subparser gets -v/-q, in one loop rather than 20 edits. Also on the top-level parser,
    # so the flag works on either side of the subcommand name. `cassette` has sub-verbs, so the loop
    # descends one level: without that, `mhcmatch cassette select -v` is an unrecognised argument
    # while `mhcmatch cassette -v select` works, which is not a distinction anybody would guess.
    _add_verbosity(ap)
    for _p in sub.choices.values():
        _add_verbosity(_p)
        for _act in _p._actions:
            if isinstance(_act, argparse._SubParsersAction):
                for _q in _act.choices.values():
                    _add_verbosity(_q)

    a = ap.parse_args(argv)
    global _V
    _V = 0 if a.quiet else (2 if a.verbose else 1)
    name = f"{a.cmd} {a.cassette_cmd}" if getattr(a, "cassette_cmd", None) else a.cmd
    if getattr(a, "deprecated", None):
        print(f"# `mhcmatch {a.cmd}` is deprecated and will be removed after 1.x; "
              f"use `mhcmatch {a.deprecated}`", file=sys.stderr)
    with step(f"mhcmatch {name}", level=2):
        a.fn(a)


if __name__ == "__main__":
    main()
