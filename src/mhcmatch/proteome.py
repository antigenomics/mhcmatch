"""Near-exact source-peptide lookup against a reference proteome.

Given a query peptide (e.g. a neoantigen), find the nearly-exact self peptide it derives from and
its parent protein / position via **full-sequence** (unmasked) ``<= max_subs`` search over all
windows of the proteome of the query's length -- using the seqtree Hamming fast path. This is a
*distinct* mode from the anchor-masked TCR-facing homology and the presentation-signature searches.
See the theory appendix §5 (near-exact source identification).
"""
from __future__ import annotations

import gzip
from array import array
from bisect import bisect_right
from dataclasses import dataclass

from seqtree import Index, SearchParams

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
    protein: str
    position: int       # 0-based start in the protein
    ref_peptide: str
    n_subs: int
    mutations: tuple    # ((pos_in_peptide, query_aa, ref_aa), ...)


class _Meta:
    """``ref_id -> (protein, position, window)`` for one window length, stored as two int arrays.

    A list of tuples is the obvious representation and it does not fit: see :meth:`Proteome._index`.
    ``starts[k]`` is the first ``ref_id`` belonging to protein ``names[k]``, so the owning protein is
    a bisect and the window is a slice of the sequence the object already holds."""

    __slots__ = ("seqs", "names", "starts", "pos", "L")

    def __init__(self, seqs, names, starts, pos, L):
        self.seqs, self.names, self.L = seqs, names, L
        self.starts = array("l", starts)
        self.pos = pos if isinstance(pos, array) else array("l", pos)

    def __len__(self):
        return len(self.pos)

    def __getitem__(self, ref_id):
        k = bisect_right(self.starts, ref_id) - 1
        name = self.names[k]
        i = self.pos[ref_id]
        return name, i, self.seqs[name][i:i + self.L].upper()


class Proteome:
    """A reference proteome with lazily-built per-length window indices."""

    def __init__(self, seqs):
        self.seqs = seqs
        self._cache = {}   # length -> (Index | None, [(protein, pos, window), ...])

    @classmethod
    def from_fasta(cls, path):
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

    def _index(self, L):
        """``(Index, meta)`` for length ``L``, built once and cached.

        ``meta`` is a ``_Meta``, not a list of ``(protein, pos, window)`` tuples, and on a whole
        proteome that is the difference between usable and not. The human proteome has **68,389,335**
        9-mer windows -- every position of every protein, not distinct sequences, because the point
        of this index is *where* a peptide comes from. One Python tuple per window holding a name
        reference, an int and a fresh 9-character string costs ~200 bytes: **~14 GB per length**, so
        a query set spanning 8-11 asked for ~55 GB of metadata alone. An ``array("l")`` of positions
        plus the protein names already in :attr:`seqs` costs 8 bytes per window and reconstructs the
        same triple on demand.

        What remains is inherent: the window list handed to ``Index.build`` and the index itself.
        Measured on the human proteome, 12.6 GB peak for the first length and ~3.6 GB for each
        further one. **Ask for the lengths you need** -- :meth:`find_sources` builds one index per
        distinct query length, so a mixed 8-11 query set builds four."""
        if L not in self._cache:
            names, starts, pos = [], [], array("l")
            windows = []
            for name, seq in self.seqs.items():
                s = seq.upper()
                names.append(name)
                starts.append(len(windows))
                for i in range(len(s) - L + 1):
                    w = s[i:i + L]
                    if all(c in _AA for c in w):
                        windows.append(w)
                        pos.append(i)
            meta = _Meta(self.seqs, names, starts, pos, L)
            self._cache[L] = (Index.build(windows, alphabet="aa") if windows else None, meta)
        return self._cache[L]

    def find_source(self, peptide, max_subs=1, exclude_exact=False):
        """Self peptides within ``max_subs`` substitutions of ``peptide``, nearest first.

        Returns ``[SourceHit, ...]``. ``exclude_exact=True`` drops perfect (0-mismatch) matches --
        useful to find the wild-type a mutated neoantigen derives from when the query is itself self.
        """
        q = peptide.strip().upper()
        idx, meta = self._index(len(q))
        if idx is None:
            return []
        p = SearchParams(max_subs=max_subs, engine="seqtm")
        out = []
        for hit in idx.search(q, p):
            name, pos, w = meta[hit.ref_id]
            muts = tuple((i, q[i], w[i]) for i in range(len(q)) if q[i] != w[i])
            if exclude_exact and not muts:
                continue
            out.append(SourceHit(name, pos, w, len(muts), muts))
        out.sort(key=lambda h: h.n_subs)
        return out

    def find_sources(self, peptides, max_subs=1, exclude_exact=False, threads=0):
        """``{peptide: [SourceHit, ...]}`` for many peptides at once -- the batch form of
        :meth:`find_source`.

        One index build per distinct length and **one threaded C++ batch query** per length
        (``search_batch`` releases the GIL), instead of one Python-level query per peptide. The index
        build dominates a single lookup -- roughly a minute for the human proteome -- so the
        per-peptide entry point is the wrong one for anything but an interactive question.

        ``threads=0`` uses every core. Duplicate and blank queries are collapsed; the returned dict
        is keyed by the stripped, upper-cased peptide."""
        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        out = {q: [] for q in qs}
        by_len = {}
        for q in qs:
            by_len.setdefault(len(q), []).append(q)
        p = SearchParams(max_subs=max_subs, engine="seqtm")
        for L, group in by_len.items():
            idx, meta = self._index(L)
            if idx is None:
                continue
            for q, hits in zip(group, idx.search_batch(group, p, threads)):
                res = []
                for hit in hits:
                    name, pos, w = meta[hit.ref_id]
                    muts = tuple((i, q[i], w[i]) for i in range(L) if q[i] != w[i])
                    if exclude_exact and not muts:
                        continue
                    res.append(SourceHit(name, pos, w, len(muts), muts))
                res.sort(key=lambda h: h.n_subs)
                out[q] = res
        return out

    def find_exact_sources(self, peptides):
        """``{peptide: [SourceHit, ...]}`` at **exactly zero substitutions** -- :meth:`find_sources`
        with ``max_subs=0``, without building the fuzzy index.

        Same return shape and the same ``(protein, position, ref_peptide, n_subs, mutations)``
        content; ``n_subs`` is 0 and ``mutations`` is ``()`` for every hit, because that is what an
        exact match is. Peptides with no source come back with an empty list, so the dict is keyed by
        every distinct stripped, upper-cased query.

        **Why a second entry point rather than a flag.** :meth:`_index` is a Python loop over every
        position of every protein -- 68,398,087 iterations for ``L = 9``, ~12.6 GB peak -- and what
        it buys is the ability to answer ``<= max_subs``. At ``max_subs = 0`` none of that is
        queried: the question is set membership, which :meth:`window_array` already answers
        vectorised in 11.0 s (its own docstring measures it), and provenance is then needed only for
        the handful of peptides that actually hit.

        Membership *and* provenance come out of one sorted array per length
        (:meth:`_sorted_windows`): every standard-AA window with the buffer offset it came from,
        sorted by sequence. A pair of ``np.searchsorted`` calls over the whole query batch gives each
        query the half-open slice of that array holding **all** of its occurrences, and one more
        ``searchsorted`` on the protein start offsets names the owning protein. Nothing scans, and
        nothing is done per residue in Python.

        **A per-hit ``bytes.find`` over the proteome buffer is the obvious alternative and it does
        not scale** -- measured, not assumed. A 27-mer unit is native context by design, so most of
        its ~70 registers *are* self peptides: on 500 units at ``L = 9``, 9,497 of 9,500 distinct
        registers hit, and a full-buffer scan apiece is O(hits x proteome) rather than
        O(hits x log windows).

        The safety screen is the caller this exists for: it runs at ``max_subs=0`` by design
        (:func:`mhcmatch.vector.self_origin_risk` measures why), over tens of thousands of registers,
        and it never asked a fuzzy question.
        """
        import numpy as np

        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        out: dict = {q: [] for q in qs}
        by_len: dict = {}
        for q in qs:
            by_len.setdefault(len(q), []).append(q)
        for L, group in by_len.items():
            win, off = self._sorted_windows(L)
            if win.size == 0:
                continue
            _buf, names, starts = self._buffer()
            q = np.array(group, dtype=f"S{L}")
            lo = np.searchsorted(win, q, side="left")
            hi = np.searchsorted(win, q, side="right")
            for k in np.flatnonzero(hi > lo):
                at = off[lo[k]:hi[k]]
                owner = np.searchsorted(starts, at, side="right") - 1
                pep = group[int(k)]
                out[pep] = [SourceHit(names[int(j)], int(a) - int(starts[j]), pep, 0, ())
                            for a, j in zip(at, owner)]
        return out

    def _sorted_windows(self, L):
        """``(windows, buffer offsets)`` for every standard-AA length-``L`` window, sorted by
        sequence and cached.

        :meth:`window_array` is the same construction with ``np.unique`` on the end, which is what a
        membership question wants and what a *provenance* question cannot use -- ``unique`` is
        exactly the information about where each window came from. The sort is stable, so the
        occurrences of one window stay in buffer order and therefore in the order
        :meth:`find_sources` reports them.

        Offsets are ``int32``: the human proteome is ~7.4e7 residues, three orders of magnitude
        below the type's range, and the array has one entry per window.
        """
        import numpy as np
        key = ("sorted", L)
        if key not in self._cache:
            buf, _names, _starts = self._buffer()
            a = np.frombuffer(buf, dtype=np.uint8)
            if a.size < L:
                self._cache[key] = (np.empty(0, dtype=f"S{L}"), np.empty(0, dtype=np.int32))
                return self._cache[key]
            ok = np.zeros(256, dtype=bool)
            for c in _AA:
                ok[ord(c)] = True
            sw = np.lib.stride_tricks.sliding_window_view
            keep = np.flatnonzero(sw(ok[a], L).all(axis=1)).astype(np.int32)
            win = np.ascontiguousarray(sw(a, L)[keep]).view(f"S{L}").ravel()
            order = np.argsort(win, kind="stable")
            self._cache[key] = (win[order], keep[order])
        return self._cache[key]

    def _buffer(self):
        """``(bytes, [protein name], np.int64 start offsets)`` -- the proteome as one NUL-joined
        buffer, cached. NUL separates proteins and is not a residue, so a window that contains only
        residues lies wholly inside one protein and ``searchsorted`` on ``starts`` names which."""
        import numpy as np
        if ("buf",) not in self._cache:
            names = list(self.seqs)
            lens = np.fromiter((len(self.seqs[n]) for n in names), dtype=np.int64,
                               count=len(names))
            starts = np.zeros(len(names), dtype=np.int64)
            if len(names):
                starts[1:] = np.cumsum(lens[:-1] + 1)
            buf = "\x00".join(self.seqs[n].upper() for n in names).encode("latin1")
            self._cache[("buf",)] = (buf, names, starts)
        return self._cache[("buf",)]

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
        """Set of all length-``L`` standard-AA proteome windows (lazy). ~1 GB/length as a Python set --
        much lighter than the seqtree index, and O(1) membership for the 1-sub wildtype fast path."""
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

        * **Only the nearest shell votes.** ``best = min(h.n_subs)``, and hits further out are
          discarded rather than pooled. A radius-2 shell is ~85x the size of the radius-1 shell
          inside it, so pooling them lets a distant coincidence outvote a genuine
          single-substitution parent.
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

        One length at a time, dropping :attr:`_cache` after each: :meth:`_index` costs 12.6 GB peak
        for the first length and ~3.6 GB for each further one, so holding a four-length query set at
        once is the difference between fitting in memory and not. ``path`` is the FASTA the symbols
        are read from and defaults to the one this proteome was loaded from, as in
        :meth:`window_genes`.
        """
        path = path or getattr(self, "_path", None)
        if path is None:
            raise ValueError("assign_genes needs the FASTA the symbols are read from")
        gene_of = gene_symbols(path, key="name")
        qs = sorted({str(p).strip().upper() for p in peptides if str(p).strip()})
        by_len = {}
        for q in qs:
            by_len.setdefault(len(q), []).append(q)
        out = {q: [] for q in qs}
        for L in sorted(by_len):
            hits = self.find_sources(by_len[L], max_subs=max_subs, exclude_exact=True,
                                     threads=threads)
            self._cache.pop(L, None)
            for q, hs in hits.items():
                if not hs:
                    continue
                best = min(h.n_subs for h in hs)
                out[q] = sorted({g for h in hs
                                 if h.n_subs == best and (g := gene_of.get(h.protein))})
        return out

    def wildtype(self, peptide, max_subs=1):
        """The wild-type self peptide a mutated ``peptide`` derives from, or ``None``.

        A self peptide exactly one substitution away (its point-mutation origin) -- the position-aligned
        WT counterpart needed for agretopicity / DAI when the caller has no WT window (e.g. a bare
        neoantigen list like TESLA). ``None`` when nothing is one sub away (indel / spliced / non-self,
        or the peptide is itself an exact self peptide with no mutated origin). Ties resolve to the
        first variant found (position, then residue order).

        For ``max_subs=1`` this uses a hash-set fast path (generate the L*19 single-sub variants and
        test proteome membership -- microseconds/peptide, so it scales to large corpora); larger
        ``max_subs`` falls back to the general :meth:`find_source` fuzzy search.
        """
        q = peptide.strip().upper()
        if max_subs == 1 and all(c in _AA for c in q):
            ws = self._window_set(len(q))
            for i in range(len(q)):
                pre, post = q[:i], q[i + 1:]
                for a in _AA_ORDER:
                    if a != q[i]:
                        v = pre + a + post
                        if v in ws:
                            return v
            return None
        hits = self.find_source(peptide, max_subs=max_subs, exclude_exact=True)
        return hits[0].ref_peptide if hits else None
