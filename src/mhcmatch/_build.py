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

#: ``(k, mask)`` the shipped tables are built at -- the face convention the scored model reads.
#: A table is a pure function of ``(cls, comp, species, k, mask)`` and nothing else, so this is the
#: whole of what a release commits to. It is a list because a version that carries two conventions
#: (a shipped one and a deprecated one still cited by a recorded result) is a real state.
#:
#: **A wildcard table is 21**k cells, not 20**k**, and the ratio bites: at ``k = 5`` that is
#: 4,084,101 float64 per table against 3,200,000, and twelve of them. Widening the shipped set is
#: therefore a size decision as well as a modelling one -- build a sweep's tables in the analysis
#: repo, and ship only the convention that won.
SHIPPED_CORPUS = [(3, "slice")]

#: Every ``(class, component, species)`` the shipped model reads. All three components are keyed on
#: species: ``self`` is a proteome and takes it explicitly, ``viral`` is one file whose own
#: ``mhc_species`` column holds both, and ``thymus`` is one file per species.
CORPUS_COMBOS = [(cls, comp, sp)
                 for cls in ("mhc1", "mhc2")
                 for comp, sp in (("thymus", "human"), ("thymus", "mouse"),
                                  ("self", "human"), ("self", "mouse"),
                                  ("viral", "human"), ("viral", "mouse"))]


def _key(cls: str, comp: str, species: str, k: int, mask: str = "slice") -> str:
    """Archive key for one table. ``slice`` keeps its 0.24-0.26 form so an old artifact still loads;
    anything else appends the mask. Must stay in step with ``mimicry._vendored_counts``."""
    stem = f"{cls}|{comp}|{species}|{k}"
    return stem if mask == "slice" else f"{stem}|{mask}"


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
    old = {}
    if os.path.exists(out):
        with np.load(out, allow_pickle=False) as z:
            old = {n: z[n] for n in z.files if n != "meta"}
    tables, meta, moved_any = {}, {}, []
    for k, mask in SHIPPED_CORPUS:
        cells = mimicry.alphabet(mask) ** k
        for cls, comp, sp in CORPUS_COMBOS:
            t0 = time.time()
            T, N = mimicry.corpus_counts(None, cls, comp, k=k, self_species=sp, mask=mask)
            dt = time.time() - t0
            name = _key(cls, comp, sp, k, mask)
            tables[name] = np.asarray(T, dtype=np.float64)
            meta[name] = {"n": float(N), "dense": int((T > 0).sum()), "seconds": round(dt, 2),
                          "k": k, "mask": mask}
            moved = ""
            if name in old and not np.array_equal(old[name], tables[name]):
                d = np.abs(old[name] - tables[name])
                moved = (f"  ** MOVED: max |old - new| = {d.max():.6g} on "
                         f"{int((d > 0).sum())} cells **")
                moved_any.append(name)
            say(f"  {cls:5s} {comp:7s} {sp:6s} k={k} {mask:8s} {dt:6.1f}s  N={N:>15,.0f}  "
                f"dense={meta[name]['dense']:>7,}/{cells:,}{moved}")
    meta["_"] = {"version": mhcmatch.__version__, "shipped": [list(x) for x in SHIPPED_CORPUS],
                 "k": SHIPPED_CORPUS[0][0], "script": "mhcmatch build corpus"}
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
                     "recognition_esm_pca.npz",
                     "recognition_default.json"]),
    # -- fitted artifacts whose generator lives in the benchmark repo -----------------
    # Until 2026-08-23 these were shipped but *unchecked*: `--check` covered 11 of the 27 files
    # in this directory, so the sixteen below could go missing or be shipped from a half-finished
    # copy and nothing would say so. They are listed here with `None` builders so `--check`
    # validates presence and `build` prints the command that regenerates them. Every command in
    # EXTERNAL below is one the artifact or PROVENANCE.md records -- none is reconstructed.
    "aggregate": ("the EPIC neoantigen scorer", None, ["aggregate_mhc1.json"]),
    "affinity": ("affinity head coefficients", None, ["affinity_mhc1.json"]),
    "potts": ("Potts affinity weights (the source of `occupancy`)", None,
              ["affinity_potts_mhc1.npz", "affinity_potts_mhc2.npz"]),
    "complement1": ("class-I complementarity heads", None,
                    ["complement_mhc1_human.json", "complement_mhc1_mouse.json"]),
    "complement2": ("class-II complementarity heads", None,
                    ["complement_mhc2_human.json", "complement_mhc2_mouse.json"]),
    "mimicry": ("the six-feature mimicry aggregate", None, ["mimicry_mhc1.json"]),
    "ligand": ("ligand flank/context span model", None, ["ligand_context.tsv"]),
    "pseudoseq": ("IMGT groove pseudosequences", None,
                  ["mhci_pseudo.fa", "mhcii_pseudo.fa"]),
    # -- static reference data: no fit, no generator recorded, presence-checked only ---
    "reference": ("static reference tables (no generator recorded)", None,
                  ["proteome_markov1.tsv", "mhc2_alpha_prior.tsv", "structure_templates.json"]),
}

#: How to build what this process cannot. Every command here is recorded by the artifact itself
#: (its ``generator`` field) or by ``PROVENANCE.md``; a target absent from this map has no
#: generator on record, and ``build`` says so rather than inventing one.
EXTERNAL = {
    "recognition": "uv run tools/build_recognition.py",
    # `aggregate_mhc1.json`'s own `generator` field. The fit writes this file directly into the
    # library checkout, so the hand-copy that let the GRAND -> EPIC rename reach the artifact but
    # not its generator is gone. `bench/run_epic.sh` runs the whole chain that leads to it.
    "aggregate": "python bench/immuno/epic_v4_fit.py --physchem rose_af5   # benchmark repo",
    "affinity": "python bench/affinity/train.py --cls mhc1 --species human   # benchmark repo",
    "potts": "python bench/affinity/fit_potts.py --cls mhc1    # and --cls mhc2; benchmark repo",
    "complement1": "python bench/neoag/complement.py --fit chowell_rebuilt --tables all",
    "complement2": "python bench/neoag/complement_mhc2.py      # benchmark repo",
    "mimicry": "python bench/neoag/mimicry_fit.py              # benchmark repo",
    "ligand": "python bench/train_spans.py                     # benchmark repo",
    "pseudoseq": "python ../tcren-ms/scripts/build_pseudo_fasta.py",
}


# ---------------------------------------------------------------- staleness

def _stamp(path: str):
    """The **package** version a shipped artifact carries, or None if it does not carry one.

    Two version vocabularies live in this directory and they are told apart **by the shape of the
    value, not by the file extension**. A *model* version is an int — EPIC is ``3``, the recognition
    heads are ``2`` — and comparing one to a package version is a category error that reports every
    head stale at every release. A *package* version is dotted (``"0.26.0"``), and an artifact that
    carries one is asserting which release built it, so it is checked.

    Blanket-exempting ``.json`` was the earlier rule, and it hid exactly one thing:
    ``mimicry_mhc1.json`` carries ``"0.12.0"`` and has gone unchecked across fifteen minor releases.

    An artifact with no version record at all is not stale — several are static reference data, and
    demanding a stamp would make the check cry wolf on files no release touches.
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
                try:
                    meta = json.loads(bytes(z["meta"]).decode())
                except (UnicodeDecodeError, json.JSONDecodeError):
                    # `meta` is not always a version record. `affinity_potts_*.npz` stores five
                    # int32 shape parameters under that name, and reading them as JSON raised --
                    # which the bare `except` below then reported as "unreadable", i.e. as a
                    # corrupt artifact. An npz that carries no version record is presence-checked
                    # only, exactly like a .json, rather than being stale at every release.
                    return None
                return meta.get("_", {}).get("version") if isinstance(meta, dict) else None
        if path.endswith(".json"):
            v = json.load(open(path)).get("version")
            # The two vocabularies, told apart by shape rather than by filename. A *model* version
            # is an int -- EPIC is 3, the recognition heads are 2 -- and comparing it to a package
            # version reports every head stale at every release. A *package* version is dotted, and
            # an artifact carrying one is making a claim about which release built it, so it is
            # checked. `mimicry_mhc1.json` is the file this distinction exists for: it carries
            # "0.12.0" and was never checked, because JSON was blanket-exempted.
            return v if isinstance(v, str) and "." in v else None
    except Exception as exc:                       # a file that will not open IS corrupt
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
