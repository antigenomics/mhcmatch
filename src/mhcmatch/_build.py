"""Release-time regeneration of the artifacts under ``mhcmatch.data``, driven by ``mhcmatch build``.

**Everything mhcmatch ships is rebuilt from the CLI, and the whole rebuild costs minutes.** There is
therefore never a reason to run a stale artifact: on a version bump, on a deposit change, or when a
definition moves, regenerate rather than reasoning about whether it matters. ``mhcmatch build
--check`` answers "is anything stale" without building anything, which is the form CI wants.

The builders live *here*, inside the package, rather than in ``tools/`` — a builder that only exists
in a source checkout is one a wheel user cannot run and, more to the point, one that gets forgotten.
``tools/build_*.py`` are kept as thin shims onto these functions.

Two of the three shipped families are pure ``mhcmatch`` and live here. The third, the recognition
heads, needs ESM2 (torch + transformers) — far heavier than the package's runtime requirements — so
it stays a PEP 723 self-contained script at ``tools/build_recognition.py`` and ``build`` reports the
exact command rather than pretending it can run it.

**A bump-only rebuild moves no prediction, and both builders are instrumented to prove rather than
assume it**: :func:`corpus_tables` prints ``** MOVED **`` for any cell that changed, and
:func:`anchor_models` is checked against the previous file by the caller. Measured for the anchor
models at 0.25.0 → 0.26.0: max |new − old| = 0 over 9,000 scorings.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import pickle
import time

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

#: Every ``(class, component, species)`` the shipped model reads. All three components are keyed on
#: species: ``self`` is a proteome and takes it explicitly, ``viral`` is one file whose own
#: ``mhc_species`` column holds both, and ``thymus`` is one file per species.
CORPUS_COMBOS = [(cls, comp, sp)
                 for cls in ("mhc1", "mhc2")
                 for comp, sp in (("thymus", "human"), ("thymus", "mouse"),
                                  ("self", "human"), ("self", "mouse"),
                                  ("viral", "human"), ("viral", "mouse"))]


def _key(cls: str, comp: str, species: str, k: int) -> str:
    return f"{cls}|{comp}|{species}|{k}"


# ---------------------------------------------------------------- targets

def anchor_models(say=print) -> list:
    """Refit and rewrite the vendored :class:`AnchorModel` pickles.

    The MHC-II register + K=3 motif EM is slow on the full corpus and a ``predict`` run triggers it
    twice, so both configs ship pre-fit and are loaded read-only — no runtime writes, so concurrent
    pipeline tasks never race on a cache.

    The load guard keys on ``mhcmatch.__version__``, so a bump without a rebuild ships models the
    library refuses to load. That is a *provenance* guard, not a correctness one: ``panel_sha`` and
    ``params`` are unchanged by a bump and the refit is deterministic.
    """
    import mhcmatch
    from .diffusion import _VENDORED_MODELS, save_vendored_anchor_model

    classes = tuple(dict.fromkeys(cls for cls, _, _ in _VENDORED_MODELS))
    store = mhcmatch.Store.from_pmhc(tier="full", species="human", classes=classes)
    written = []
    for (cls, footprint, background), name in _VENDORED_MODELS.items():
        path = os.path.join(DATA, name)
        t0 = time.time()
        save_vendored_anchor_model(store, cls, path, footprint=footprint, background=background)
        written.append(path)
        say(f"  {name}  {os.path.getsize(path) / 1e6:.2f} MB  "
            f"[{cls} {footprint}/{background}]  {time.time() - t0:.1f} s")
    return written


def corpus_tables(say=print) -> list:
    """Rebuild the vendored corpus k-mer count tables.

    ``mimicry.corpus_counts`` slides over every window of a reference set; for ``self`` that is the
    whole proteome — ~12.7 M windows per length, four lengths at class I — measured at 53.0 s (mhc1)
    and 14.0 s (mhc2) *per process*, before a single peptide is scored. The result is 8,000 float64s.
    So it ships.

    Unlike the anchor models, a bump-only rebuild here is genuinely bit-identical, and that is
    asserted rather than assumed: any cell that moves is printed as ``** MOVED **`` at build time
    instead of surfacing later as a shifted AUROC.
    """
    import numpy as np

    import mhcmatch
    from . import mimicry

    out = os.path.join(DATA, "corpus_tables.npz")
    # **Build, do not read.** `corpus_counts` serves the vendored artifact on the default path, so a
    # builder that did not disable it would re-emit whatever is already committed with a fresh
    # version stamp -- the MOVED check could never fire and a deposit change would ship silently.
    # Empty dict, not None: None means "not loaded yet" and would fall back to disk.
    mimicry._VENDORED = {}
    mimicry._COUNTS.clear()
    k = mimicry.CORPUS_K
    old = {}
    if os.path.exists(out):
        with np.load(out, allow_pickle=False) as z:
            old = {n: z[n] for n in z.files if n != "meta"}
    tables, meta, moved_any = {}, {}, []
    for cls, comp, sp in CORPUS_COMBOS:
        t0 = time.time()
        T, N = mimicry.corpus_counts(None, cls, comp, k=k, self_species=sp)
        dt = time.time() - t0
        name = _key(cls, comp, sp, k)
        tables[name] = np.asarray(T, dtype=np.float64)
        meta[name] = {"n": float(N), "dense": int((T > 0).sum()), "seconds": round(dt, 2)}
        moved = ""
        if name in old and not np.array_equal(old[name], tables[name]):
            d = np.abs(old[name] - tables[name])
            moved = f"  ** MOVED: max |old - new| = {d.max():.6g} on {int((d > 0).sum())} cells **"
            moved_any.append(name)
        say(f"  {cls:5s} {comp:7s} {sp:6s} {dt:6.1f}s  N={N:>15,.0f}  "
            f"dense={meta[name]['dense']:>5,}/{20 ** k:,}{moved}")
    meta["_"] = {"version": mhcmatch.__version__, "k": k, "script": "mhcmatch build corpus"}
    buf = io.BytesIO()
    np.savez_compressed(buf, meta=np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
                        **tables)
    with open(out, "wb") as fh:
        fh.write(buf.getvalue())
    say(f"  wrote {os.path.basename(out)}  ({os.path.getsize(out) / 1e3:.1f} kB, "
        f"{len(tables)} tables)")
    if moved_any:
        say(f"  ** {len(moved_any)} table(s) moved: {', '.join(moved_any)} -- this is a data change, "
            f"not a version bump; re-baseline anything downstream **")
    return [out]


#: ``target -> (human name, builder or None, the files it owns)``. ``None`` means the artifact is
#: real and shipped but cannot be built in-process; ``--check`` still validates it and ``build``
#: prints the command that does build it.
TARGETS = {
    "anchor": ("vendored AnchorModel pickles", anchor_models,
               ["anchor_model_mhc1_proteome_adaptive.pkl.gz",
                "anchor_model_mhc2_proteome_adaptive.pkl.gz",
                "anchor_model_mhc2_proteome_core.pkl.gz"]),
    "corpus": ("corpus k-mer count tables", corpus_tables, ["corpus_tables.npz"]),
    "recognition": ("recognition heads (needs ESM2)", None,
                    ["recognition_esm64_glm_mhc1_human.json",
                     "recognition_esm64_glm_mhc1_mouse.json",
                     "recognition_physchem_glm_mhc1_human.json",
                     "recognition_physchem_glm_mhc1_mouse.json",
                     "recognition_posbayes_mhc1_human.json",
                     "recognition_posbayes_mhc1_mouse.json",
                     "recognition_esm_pca.npz"]),
}

#: How to build what this process cannot.
EXTERNAL = {"recognition": "uv run tools/build_recognition.py"}


# ---------------------------------------------------------------- staleness

def _stamp(path: str):
    """The **package** version a shipped artifact carries, or None if it does not carry one.

    Only the pickles and the npz stamp ``mhcmatch.__version__``; those are the ones whose load guards
    key on it. A ``.json`` artifact's ``version`` is a *model* version — EPIC is ``3``, the
    recognition heads are ``2`` — and comparing that to a package version is a category error that
    reports every head stale at every release. So JSON returns None here and is checked for presence
    only.
    """
    import numpy as np
    try:
        if path.endswith(".pkl.gz"):
            meta, _ = pickle.loads(gzip.decompress(open(path, "rb").read()))
            return meta.get("version")
        if path.endswith(".npz"):
            with np.load(path, allow_pickle=False) as z:
                if "meta" not in z.files:
                    return None
                return json.loads(bytes(z["meta"]).decode()).get("_", {}).get("version")
    except Exception as exc:                       # a corrupt artifact is stale by definition
        return f"unreadable: {type(exc).__name__}"
    return None


def check(targets=None) -> list:
    """Every shipped artifact's version stamp against ``__version__``.

    Returns the rows that are **stale** — each ``(target, filename, found, want)``. A stamp of
    ``None`` is not stale: several artifacts are static reference data that carries no version, and
    demanding one would make the check cry wolf on files no release touches.
    """
    from . import __version__
    bad = []
    for name in (targets or TARGETS):
        _, _, files = TARGETS[name]
        for f in files:
            p = os.path.join(DATA, f)
            if not os.path.exists(p):
                bad.append((name, f, "MISSING", __version__))
                continue
            got = _stamp(p)
            if got is not None and got != __version__:
                bad.append((name, f, got, __version__))
    return bad
