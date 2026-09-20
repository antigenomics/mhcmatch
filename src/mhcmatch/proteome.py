"""Near-exact source-peptide lookup against a reference proteome.

Given a query peptide (e.g. a neoantigen), find the nearly-exact self peptide it derives from and
its parent protein / position via **full-sequence** (unmasked) ``<= max_subs`` search over the
proteome -- using :class:`seqtree.TextIndex`. This is a *distinct* mode from the anchor-masked
TCR-facing homology and the presentation-signature searches. See the theory appendix §5 (near-exact
source identification).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass

from seqtree import TextIndex

#: Seed width of the one index this module builds. ``k`` belongs to the **index**, not to the query,
#: which is the whole reason nothing here is per length any more: one build answers every peptide
#: length and every radius.
#:
#: 4 rather than 5, for two reasons. A query shorter than ``k`` cannot be answered at all, and the
#: shortest register we ask about is an 8-mer. And the search splits a query into
#: ``min(max_subs + 1, L // k)`` blocks, so at radius 2 a 9-mer gets two blocks at ``k = 4`` and only
#: one at ``k = 5`` -- fewer blocks filter worse. Measured on the human proteome (147,506 records,
#: 69,578,135 residues): 0.7 s build / 0.6 GB peak RSS at ``k = 4``, 0.7 s / 1.0 GB at ``k = 5``.
_TEXT_K = 4

_AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_AA = set(_AA_ORDER)


def read_fasta(path):
    """``{name: sequence}`` from a (optionally gzipped) FASTA; name = first whitespace token."""
    op = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def gene_symbols(path, key: str = "name"):
    """``{key: gene}`` from the UniProt ``GN=`` field. ``key="name"`` (default) matches
    :func:`read_fasta`; ``key="accession"`` matches a bare UniProt accession.

    **Both keyings exist because two different callers need two different sides of the same header.**
    A :class:`SourceHit` names its protein as the FASTA's first whitespace token,
    ``sp|Q8WZ42|TITIN_HUMAN``, so a proteome scan needs ``name``. The thymic and ligandome deposits
    record ``source_protein`` as a bare accession, ``Q8WZ42``, so
    :func:`mhcmatch.mimicry.safety` needs ``accession``. Neither can reach
    :func:`mhcmatch.expression.safety_profile`, which is keyed on the HGNC symbol ``TTN``, without
    one of them.

    Without it there is no way to ask *which tissue* a T cell cross-reactive with a given self peptide
    would attack, and that is the question separating a titin match (``Q8WZ42`` → ``TTN`` → heart left
    ventricle, 64 TPM) from a testis-restricted one.

    The symbol is absent from :func:`read_fasta`'s output because that function keeps only the name,
    and widening its return contract would ripple through every caller. A second pass over the
    headers is cheap -- one second for the human proteome -- and additive.

    Entries with no ``GN=`` map to ``None`` rather than being dropped: the 147,506 human records
    include TrEMBL entries with no assigned symbol, and silently losing them would overstate the
    coverage of any downstream tissue filter.
    """
    import re

    if key not in ("name", "accession"):
        raise ValueError(f"key must be 'name' or 'accession', got {key!r}")
    op = gzip.open if str(path).endswith(".gz") else open
    pat = re.compile(r"\bGN=(\S+)")
    out = {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].rstrip()
                name = head.split()[0]
                m = pat.search(head)
                gene = m.group(1) if m else None
                if key == "name":
                    out[name] = gene
                else:
                    # sp|Q8WZ42|TITIN_HUMAN -> Q8WZ42; a header without the db|acc|id form is
                    # keyed on itself rather than dropped.
                    parts = name.split("|")
                    out[parts[1] if len(parts) >= 3 else name] = gene
    return out


@dataclass
class SourceHit:
    """One near-exact match of a query peptide against a reference window -- the result row of
    :meth:`Proteome.find_source` / :meth:`find_sources` / :meth:`find_exact_sources`."""
    protein: str
    position: int       # 0-based start in the protein
    ref_peptide: str
    n_subs: int
    mutations: tuple    # ((pos_in_peptide, query_aa, ref_aa), ...)


class Proteome:
    """A reference proteome, searched by peptide through one :class:`seqtree.TextIndex`."""

    def __init__(self, seqs):
        self.seqs = seqs
        #: ``ref_id -> name``. `TextIndex` addresses records positionally -- by their index in the
        #: list handed to ``build`` -- so this is that list, and its order is load-bearing.
        self._names = list(seqs)
        self._text = None
        self._cache = {}   # memo for the window views below, keyed by (kind, L)

    @classmethod
    def from_fasta(cls, path):
        """A :class:`Proteome` over a caller-supplied FASTA, overriding :meth:`from_hf`'s fetch."""
        pm = cls(read_fasta(path))
        pm._path = path            # so `window_genes` can re-read the headers for GN=
        return pm

    @classmethod
    def from_hf(cls, name="human"):
        """Load a reference proteome by name, auto-fetched from the public HF dataset (no manual
        download). ``name`` = ``"human"`` / ``"mouse"`` (UP000005640 / UP000000589) or a pathogen
        stem; see :func:`mhcmatch.store.fetch_proteome`."""
        from .store import fetch_proteome
        return cls.from_fasta(fetch_proteome(name))

    def _index(self):
        """The one :class:`seqtree.TextIndex` over this proteome, built once per instance.

        **It takes no length, and that is the entire point.** The index this replaced was built per
        query length over every window of every protein -- 68,389,335 of them at ``L = 9`` for the
        human proteome, ~65 s and 12.6 GB peak for the first length and ~3.6 GB for each further
        one, so a class-II query set spanning 11-25 asked for fifteen of them (~82 GB). `TextIndex`
        indexes the sequences themselves and puts ``k`` in the index rather than in the query, so one
        build answers every length and every radius: measured on the same human proteome,
        **0.7 s and 0.6 GB peak, once**.

        That is also why nothing caches this to disk and nothing hands off between concurrent
        builders. Both existed only to amortise the old build -- a ``blake2b`` key over every
        sequence, a ``.building`` marker claimed with ``O_EXCL``, a bounded wait with a backing-off
        poll, and a cache directory nothing ever evicted, measured at 225 GB. **The race-free design
        is the one with nothing to race on**, and at 0.7 s rebuilding is cheaper than agreeing about
        who rebuilds.
        """
        if self._text is None:
            self._text = TextIndex.build([s.upper() for s in self.seqs.values()],
                                         alphabet="aa", k=_TEXT_K)
        return self._text

    def _rows(self, q, hits):
        """`TextIndex` hits for one query as `SourceHit`s, dropping any non-standard window.

        ``h.offset`` is already the 0-based start within the record, and ``h.mismatches`` is already
        ``(position_in_query, query_aa, ref_aa)`` in ascending position -- exactly what
        :class:`SourceHit` carries, and exactly what the ``tuple((i, q[i], w[i]) ...)`` comprehension
        this replaced produced. The window is sliced out of :attr:`seqs` rather than out of
        ``TextIndex.ref_seq``, which is lossy: it reports any residue outside its codec as ``X``.

        **The standard-AA filter is not cosmetic.** The index this replaced dropped any window
        holding a non-standard residue outright, where `TextIndex` keeps ``X`` as a real, searchable
        symbol -- only ``U`` (36 residues of the human proteome's 69,578,135) is outside the codec
        and becomes a hole no hit crosses. Re-applying the filter is what makes this port
        hit-for-hit identical to what it replaced.
        """
        L, rows = len(q), []
        for h in hits:
            name = self._names[h.ref_id]
            w = self.seqs[name][h.offset:h.offset + L].upper()
            if not _AA.issuperset(w):
                continue
            rows.append(SourceHit(name, h.offset, w, h.n_subs, tuple(h.mismatches)))
        return rows

    def _hits(self, queries, max_subs, exclude_exact=False, best_only=False, threads=0):
        """``{query: [SourceHit, ...]}``, nearest first -- the one path from a hit to a `SourceHit`.

        **``best_only`` and the standard-AA filter do not compose, and getting that wrong is
        silent.** The search picks the nearest non-empty shell *before* :meth:`_rows` sees it, so a
        shell made entirely of windows holding a non-standard residue comes back empty rather than
        falling through to the next shell -- and the query then looks like it has no parent at all.
        Measured on the real corpus: **1 peptide of 445,466**, the 8-mer ``VHTCSPTN``, whose nearest
        parent is a single ``X``-bearing window at one substitution, with eight clean parents behind
        it at two. It resolved to five tied genes before this and to nothing after, which is a
        *worse* answer than the slow one and the exact failure this port had to avoid.

        So a query the filter empties is re-asked without ``best_only`` and reduced here. It is a
        second pass over a handful of queries, not over the batch.
        """
        out: dict = {q: [] for q in queries}
        qs = [q for q in queries if _AA.issuperset(q)]
        if not qs:
            return out
        ix = self._index()
        res = ix.search_batch(qs, max_subs=max_subs, exclude_exact=exclude_exact,
                              best_only=best_only, threads=threads)
        retry = []
        for q, hits in zip(qs, res):
            rows = self._rows(q, hits)
            if best_only and not rows and len(hits):
                retry.append(q)
                continue
            rows.sort(key=lambda r: r.n_subs)
            out[q] = rows
        if retry:
            res = ix.search_batch(retry, max_subs=max_subs, exclude_exact=exclude_exact,
                                  threads=threads)
            for q, hits in zip(retry, res):
                rows = self._rows(q, hits)
                if rows:
                    best = min(r.n_subs for r in rows)
                    rows = [r for r in rows if r.n_subs == best]
                out[q] = rows
        return out

    def find_source(self, peptide, max_subs=1, exclude_exact=False):
        """Self peptides within ``max_subs`` substitutions of ``peptide``, nearest first.

        Returns ``[SourceHit, ...]``. ``exclude_exact=True`` drops perfect (0-mismatch) matches --
        useful to find the wild-type a mutated neoantigen derives from when the query is itself self.

        The single-query form of :meth:`find_sources`. It used to be the wrong entry point for
        anything but an interactive question, because the index build dominated it; the index is now
        one shared build of 0.7 s, so asking about one peptide costs one query.
        """
        q = peptide.strip().upper()
        return self._hits([q], max_subs, exclude_exact)[q]

    def find_sources(self, peptides, max_subs=1, exclude_exact=False, threads=0, best_only=False):
        """``{peptide: [SourceHit, ...]}`` for many peptides at once -- the batch form of
        :meth:`find_source`.

        **One index and one threaded C++ batch query for the whole set**, whatever lengths it spans:
        ``search_batch`` releases the GIL, and ``threads=0`` uses every core. The per-length loop
        this replaced is gone, and with it the advice to ask only for the lengths you need -- there
        is nothing left here that is per length.

        ``best_only=True`` returns only the nearest non-empty shell. That is what a caller voting on
        a parent wants (:meth:`assign_genes`), and asking the search for it is far cheaper than
        materialising the outer shells in order to discard them: measured on 487,000 peptides at
        radius 2, **2.2 s for 6,824,481 hits against 21.6 s for 15,822,792**.

        Duplicate and blank queries are collapsed; the returned dict is keyed by the stripped,
        upper-cased peptide. A query carrying a residue outside the 20 standard amino acids maps to
        ``[]``, because every window that could have matched one is excluded from the answer anyway.
        """
        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        return self._hits(qs, max_subs, exclude_exact, best_only, threads)

    def find_exact_sources(self, peptides):
        """``{peptide: [SourceHit, ...]}`` at **exactly zero substitutions** -- :meth:`find_sources`
        with ``max_subs=0``.

        Same return shape and the same ``(protein, position, ref_peptide, n_subs, mutations)``
        content; ``n_subs`` is 0 and ``mutations`` is ``()`` for every hit, because that is what an
        exact match is. Peptides with no source come back with an empty list, so the dict is keyed by
        every distinct stripped, upper-cased query.

        **It survives as a name, not as a separate implementation.** It existed because the index it
        would otherwise have queried was a Python loop over every position of every protein (~12.6 GB
        peak at ``L = 9``), buying an ability an exact question never uses -- so it answered
        membership out of a sorted window array and a pair of ``np.searchsorted`` calls instead, in
        11.0 s. One `TextIndex` build is 0.7 s and answers radius 0 as directly as any other radius,
        so the second construction bought nothing and is gone. The name stays because
        :func:`mhcmatch.vector.self_origin_risk` dispatches on it, and that dispatch is how the code
        says the safety screen asks an exact question.
        """
        return self.find_sources(peptides, max_subs=0)

    def windows(self, L):
        """Every distinct length-``L`` standard-AA window of the proteome, as a ``set``. Cached, and
        roughly 1 GB per length for the human proteome -- ask for the lengths you need."""
        return self._window_set(L)

    def window_array(self, L):
        """Every distinct length-``L`` standard-AA window, as a **sorted** ``|S{L}`` numpy array.

        The vectorized form of :meth:`windows`, and 2.7x faster on the human proteome: 11.0 s
        against 30.0 s for 12,073,995 distinct 9-mers, identical output. The loop it replaces ran
        ``all(c in _AA for c in w)`` per window -- 12 M windows x 9 residues of Python-level
        membership tests.

        The array form is what a consumer that projects or indexes these wants; :meth:`windows`
        still returns a ``set`` for the O(1) membership its own callers need, and materialising a
        12 M-element Python set is most of the cost that buys.

        Packing the residues into ``uint64`` (5 bits each, ``L <= 12``) and sorting integers instead
        was tried and is **4x slower** -- 44.5 s -- because the shift/or loop costs more than numpy's
        fixed-width byte sort saves. Measured, not assumed.
        """
        import numpy as np
        key = ("arr", L)
        if key not in self._cache:
            ok = np.zeros(256, dtype=bool)
            for c in _AA:
                ok[ord(c)] = True
            # one contiguous buffer; NUL between proteins so no window straddles two, and NUL is
            # not a standard residue so the same mask rejects both cases in one pass.
            buf = np.frombuffer("\x00".join(v.upper() for v in self.seqs.values()).encode("latin1"),
                                dtype=np.uint8)
            if buf.size < L:
                self._cache[key] = np.empty(0, dtype=f"S{L}")
                return self._cache[key]
            sw = np.lib.stride_tricks.sliding_window_view
            keep = np.flatnonzero(sw(ok[buf], L).all(axis=1))
            V = np.ascontiguousarray(sw(buf, L)[keep])
            self._cache[key] = np.unique(V.view(f"S{L}").ravel())
        return self._cache[key]

    def _window_set(self, L):
        """Set of all length-``L`` standard-AA proteome windows (lazy). ~1 GB/length as a Python
        set, for callers that need O(1) membership rather than provenance.

        **Membership, not origin.** :meth:`wildtype` used to answer out of this, by generating a
        peptide's ``L * 19`` single-substitution variants and testing each -- it asks
        :meth:`find_sources` now, which needs no per-length structure at all."""
        key = ("set", L)
        if key not in self._cache:
            s = set()
            for seq in self.seqs.values():
                seq = seq.upper()
                for i in range(len(seq) - L + 1):
                    w = seq[i:i + L]
                    if all(c in _AA for c in w):
                        s.add(w)
            self._cache[key] = s
        return self._cache[key]

    def window_genes(self, peptides, path=None):
        """``{peptide: gene_symbol}`` for those of ``peptides`` that are proteome windows.

        The question a neoantigen table has to answer before it can look up expression: a candidate
        carries a somatic substitution, so it is not itself a window, but **its wild type is** --
        and that window names the gene. :meth:`wildtype` supplies the germline counterpart and this
        supplies the symbol, which is what GTEx and TCGA are keyed on.

        Streams :attr:`seqs` once and keeps only matching windows, rather than indexing the whole
        proteome and querying it: the query set is known in advance and small (~350k) where the
        index is ~68M windows per length. ``path`` is the FASTA the symbols are read from with
        :func:`gene_symbols`; it defaults to the one this proteome was loaded from.
        """
        path = path or getattr(self, "_path", None)
        if path is None:
            raise ValueError("window_genes needs the FASTA the symbols are read from")
        gene_of = gene_symbols(path, key="name")
        keys, lens = set(peptides), sorted({len(p) for p in peptides})
        out = {}
        for name, seq in self.seqs.items():
            g = gene_of.get(name)
            if not g:
                continue
            seq = seq.upper()
            for L in lens:
                for i in range(len(seq) - L + 1):
                    w = seq[i:i + L]
                    if w in keys:
                        out.setdefault(w, g)
        return out

    def assign_genes(self, peptides, max_subs=2, threads=0, path=None):
        """``{peptide: [gene, ...]}`` -- the HGNC symbol(s) of the gene each peptide derives from.

        :meth:`window_genes` answers this for a peptide that *is* a proteome window. A neoantigen is
        not: it carries the substitution that made it one, so it has to be found by near-exact
        search. That is what makes this the entry point a candidate table needs --
        ``expression.gene_level`` and both fitted expression terms are keyed on the symbol, and a
        row without one contributes a single mean-imputed constant. Measured over the benchmark
        corpus, **356,387 of 695,811 rows (51.2%) and 5,205 of 5,833 positives (89.2%)** carried no
        deposited symbol; on the VACCIMEL screen ``expr_norm`` had standard deviation **exactly
        0.0000** and AUROC **exactly 0.5000** while carrying the largest positive coefficient of the
        then-shipped EPIC v10 artifact, **+0.4950** log-odds per standard deviation. Repairing
        the symbol is what took that term to **+0.2155** in v11, on a measurement rather than a
        constant.

        Three choices, all load-bearing:

        * **Only the nearest shell votes.** A radius-2 shell is ~85x the size of the radius-1 shell
          inside it, so pooling them lets a distant coincidence outvote a genuine single-substitution
          parent. This is ``best_only=True`` -- a property of the search, not a filter over its
          output, and 10x cheaper than the filter it used to be (:meth:`find_sources`).
        * **Exact matches are excluded.** A peptide that *is* a proteome window is not a
          neoantigen, and its own gene is not the question being asked.
        * **Ties come back in full, sorted.** Resolving one needs expression data this method does
          not have -- the caller picks among the tied genes (the CLI emits a row each and lets the
          scorer's best-per-peptide selection decide). A peptide with no parent, or whose parents
          carry no ``GN=``, maps to ``[]``; neither is an error.

        ``max_subs`` defaults to 2 because **a neoantigen can carry more than one mutation**, and
        the radius is what buys the coverage: at radius 1 VACCIMEL resolves 88.2% of its peptides,
        at radius 2 96.8% (TESLA 98.5%, ITSNdb 99.5%, GBM 94.0%, Sahin_TNBC 100%).
        ``bench/results/gene_resolution.md``.

        **Nothing here is per hit, and at this scale that is the difference that matters.** The
        nearest shell of a 487,000-peptide radius-2 query is 6,824,481 hits; one `SourceHit` dataclass
        apiece is ~1.4 GB of Python objects for a function that only ever reads ``(n_subs, protein)``.
        So this reads the flat arrays out of ``search_batch`` directly and maps them through a
        ``ref_id -> gene`` integer table built once over the proteome. ``path`` is the FASTA the
        symbols are read from, and defaults to the one this proteome was loaded from, as in
        :meth:`window_genes`.
        """
        path = path or getattr(self, "_path", None)
        if path is None:
            raise ValueError("assign_genes needs the FASTA the symbols are read from")
        import numpy as np

        gene_of = gene_symbols(path, key="name")
        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        out: dict = {q: [] for q in qs}
        ok = [q for q in qs if _AA.issuperset(q)]
        if not ok:
            return out
        res = self._index().search_batch(ok, max_subs=max_subs, exclude_exact=True,
                                         best_only=True, threads=threads)
        a = res.to_numpy()
        genes: list = []
        seen: dict = {}
        gid = np.empty(len(self._names), dtype=np.int32)
        for i, name in enumerate(self._names):
            g = gene_of.get(name)
            if g is None:
                gid[i] = -1
                continue
            j = seen.get(g)
            if j is None:
                j = seen[g] = len(genes)
                genes.append(g)
            gid[i] = j
        hit_gene = gid[a["ref_id"]]
        # **The standard-AA window filter, without slicing a single window.** A non-standard residue
        # in a matched window can only sit at a MISMATCH position -- the query is all standard, so a
        # text ``X`` never matches one -- which makes "the window is clean" exactly "no mismatch
        # reads a non-standard residue". `np.cumsum` over the mismatch CSR answers that for every
        # hit at once; `np.add.reduceat` would not, because it returns the element rather than 0 for
        # an empty segment, and an exact hit has no mismatches.
        mmb, txt = a["mm_begin"], a["mm_text_aa"]
        clean = np.ones(len(hit_gene), dtype=bool)
        if txt.size:
            std = np.frombuffer(_AA_ORDER.encode("ascii"), dtype="S1")
            c = np.concatenate(([0], np.cumsum(~np.isin(txt, std))))
            clean = c[mmb[1:]] - c[mmb[:-1]] == 0
            hit_gene = np.where(clean, hit_gene, -1)
        beg = a["query_begin"]
        # **`best_only` picks the shell before the window filter runs**, so a nearest shell made
        # entirely of non-standard windows empties to nothing rather than falling through to the
        # next one -- and the peptide then looks parentless. `_hits` re-asks those without
        # `best_only`; measured on the real corpus it is 1 peptide of 445,466, and skipping the
        # re-ask cost that peptide its five tied genes.
        retry = []
        for k, q in enumerate(ok):
            lo, hi = beg[k], beg[k + 1]
            if lo == hi:
                continue
            seg = hit_gene[lo:hi]
            if not clean[lo:hi].any():
                retry.append(q)
                continue
            seg = seg[seg >= 0]
            if seg.size:
                out[q] = sorted({genes[int(j)] for j in np.unique(seg)})
        for q, hs in self._hits(retry, max_subs, exclude_exact=True).items() if retry else ():
            out[q] = sorted({g for h in hs if (g := gene_of.get(h.protein))})
        return out

    def wildtype(self, peptide, max_subs=1):
        """The wild-type self peptide a mutated ``peptide`` derives from, or ``None``.

        A self peptide exactly one substitution away (its point-mutation origin) -- the position-aligned
        WT counterpart needed for agretopicity / DAI when the caller has no WT window (e.g. a bare
        neoantigen list like TESLA). ``None`` when nothing is one sub away (indel / spliced / non-self,
        or the peptide is itself an exact self peptide with no mutated origin). Ties resolve to the
        first variant found (position, then residue order).

        One peptide at a time. :meth:`wildtypes` is the batch form and is what a corpus should call:
        the index is shared either way, but one ``search_batch`` releases the GIL and uses every core
        where a Python loop over the same peptides does neither.
        """
        q = peptide.strip().upper()
        return self.wildtypes([q], max_subs=max_subs).get(q)

    def wildtypes(self, peptides, max_subs=1, threads=0):
        """``{peptide: wild-type | None}`` -- the batch form of :meth:`wildtype`.

        One threaded ``search_batch`` for the whole corpus, keyed by the stripped, upper-cased
        peptide. It replaces two constructions at once: a per-peptide Python loop over the caller's
        list, and the ``L * 19``-variant hash-set membership test :meth:`wildtype` used to run
        against a ``_window_set`` costing ~1 GB per length.

        **The tie order is part of the contract, and it does not survive the port for free.** The
        hash-set path walked the peptide left to right and :data:`_AA_ORDER` within each position, so
        a peptide with two one-substitution parents got the earlier position, and the earlier residue
        at that position. `TextIndex` orders hits by ``(n_subs, ref_id, offset)`` -- a different
        answer on the same data, feeding ``agretopicity`` in a shipped fit. So the radius-1 shell is
        re-sorted here rather than inherited. Beyond radius 1 there was never a documented order and
        there is none now: the nearest shell wins, and ties within it go to the index's own order.
        """
        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        hits = self.find_sources(qs, max_subs=max_subs, exclude_exact=True,
                                 best_only=True, threads=threads)
        out: dict = {}
        for q in qs:
            hs = hits[q]
            if not hs:
                out[q] = None
                continue
            best = min(h.n_subs for h in hs)
            near = [h for h in hs if h.n_subs == best]
            if best == 1:
                near.sort(key=lambda h: (h.mutations[0][0], _AA_ORDER.index(h.mutations[0][2])))
            out[q] = near[0].ref_peptide
        return out
