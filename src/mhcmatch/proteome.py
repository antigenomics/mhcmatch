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
        return cls(read_fasta(path))

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

    def windows(self, L):
        """Every distinct length-``L`` standard-AA window of the proteome, as a ``set``. Cached, and
        roughly 1 GB per length for the human proteome -- ask for the lengths you need."""
        return self._window_set(L)

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
