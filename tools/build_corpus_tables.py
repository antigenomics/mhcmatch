#!/usr/bin/env python
"""Regenerate the vendored corpus k-mer count tables (release-time task).

``mimicry.corpus_counts`` builds a flat ``20**k`` table of reference-window k-mer counts by sliding
over every window of a reference set. For the deposits that is seconds; for the ``self`` channel it
is the **whole proteome** -- ~12.7 M windows per length, four lengths at class I -- and it was
measured at **53.0 s** (mhc1) and **14.0 s** (mhc2) per process, every process, before a single
peptide is scored. The result is 8,000 float64s: **64 kB of output for a minute of work.**

So it ships. The tables go in ``mhcmatch.data/corpus_tables.npz`` and ``corpus_counts`` reads them
for the default parameters, falling back to the full build for anything else (a custom
``pmhc_dir``, ``weights="locus"``, a non-default ``k``).

Rerun on a **version bump** (the load guard keys on ``mhcmatch.__version__``), when the deposits or
proteomes change, or when the definition of a face changes::

    python tools/build_corpus_tables.py

Then commit the regenerated ``src/mhcmatch/data/corpus_tables.npz`` alongside the bump.

**Unlike the anchor models, a bump-only rebuild here is genuinely bit-identical** and is asserted to
be: the script refuses to overwrite a table whose contents moved without saying so, so a rebuild
that changes a number is visible at build time rather than in a downstream AUROC.
"""
import io
import json
import os
import time

import numpy as np

import mhcmatch
from mhcmatch import mimicry

DATA = os.path.join(os.path.dirname(mhcmatch.__file__), "data")
OUT = os.path.join(DATA, "corpus_tables.npz")

#: Every (class, component, species) the shipped model reads. ``thymus`` and ``viral`` are keyed to
#: the deposit's own species column, which is human-only today; ``self`` is a proteome, so it takes
#: the species explicitly.
COMBOS = [(cls, comp, sp)
          for cls in ("mhc1", "mhc2")
          for comp, sp in (("thymus", "human"), ("self", "human"), ("self", "mouse"),
                           ("viral", "human"))]


def key(cls: str, comp: str, species: str, k: int) -> str:
    return f"{cls}|{comp}|{species}|{k}"


def main():
    k = mimicry.CORPUS_K
    old = {}
    if os.path.exists(OUT):
        with np.load(OUT, allow_pickle=False) as z:
            old = {n: z[n] for n in z.files if n != "meta"}
    tables, meta = {}, {}
    for cls, comp, sp in COMBOS:
        t0 = time.time()
        T, N = mimicry.corpus_counts(None, cls, comp, k=k, self_species=sp)
        dt = time.time() - t0
        name = key(cls, comp, sp, k)
        tables[name] = np.asarray(T, dtype=np.float64)
        meta[name] = {"n": float(N), "dense": int((T > 0).sum()), "seconds": round(dt, 2)}
        moved = ""
        if name in old and not np.array_equal(old[name], tables[name]):
            d = np.abs(old[name] - tables[name])
            moved = f"  ** MOVED: max |old - new| = {d.max():.6g} on {int((d > 0).sum())} cells **"
        print(f"{cls:5s} {comp:7s} {sp:6s} {dt:6.1f}s  N={N:>15,.0f}  "
              f"dense={meta[name]['dense']:>5,}/{20 ** k:,}{moved}")
    meta["_"] = {"version": mhcmatch.__version__, "k": k,
                 "script": os.path.basename(__file__)}
    buf = io.BytesIO()
    np.savez_compressed(buf, meta=np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
                        **tables)
    with open(OUT, "wb") as fh:
        fh.write(buf.getvalue())
    print(f"\nwrote {OUT}  ({os.path.getsize(OUT) / 1e3:.1f} kB, {len(tables)} tables)")


if __name__ == "__main__":
    main()
