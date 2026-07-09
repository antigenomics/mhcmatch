#!/usr/bin/env python3
"""Test-set construction for the head-to-head: canonical loading, eval-allele selection, and the
two splits (per-pMHC rare/frequent holdout, and zero-shot leave-one-allele-out).

Exclusion is per-pMHC (benchmark-only), matching ``bench_diffusion``/``tune_diffusion``: only the
held (epitope, allele) pair is dropped from mhcmatch's training. NetMHCpan needs no training set --
it is pretrained -- so a rare/zero-shot allele is a *handicap against mhcmatch* (NetMHCpan very
likely trained on that allele's ligands). A mhcmatch win there is therefore strong; parity is the
expected floor. See the plan's "reality check".
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # sibling compare
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # bench/

import alleles as al  # noqa: E402
from bench_diffusion import load  # noqa: E402

from mhcmatch.pseudoseq import load_pseudo, normalize_allele  # noqa: E402

_LABEL = {"mhc1": "MHCI", "mhc2": "MHCII"}
_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus"}


def load_canonical(pmhc_dir: str, cls: str, species: str = "human", tier: str = "full") -> dict:
    """``bench_diffusion.load`` re-keyed to canonical pseudoseq keys.

    Class-I keys lose the ``*`` (``HLA-A*02:01`` -> ``HLA-A02:01``) so one key space flows to both
    mhcmatch (pseudoseq lookup) and the NetMHCpan emitter; class-II keys are already canonical."""
    path = os.path.join(pmhc_dir, "pmhc", f"pmhc_{tier}.tsv.gz")
    raw = load(path, cls, _SPECIES[species])
    if cls == "mhc2":
        return raw
    out: dict = {}
    for a, peps in raw.items():
        out.setdefault(normalize_allele(a), {}).update(peps)
    return out


def select_eval_alleles(refcount: dict, cls: str, rng, n_sample: int = 20,
                        rare_max: int = 30, freq_min: int = 200) -> list:
    """Pseudo-matched, tool-supported eval alleles: all rare + a sample of medium + of frequent.

    Returns a **sorted** list -- deterministic iteration order matters because downstream splitting
    draws the seeded RNG per allele, and set iteration over allele strings is hash-randomized."""
    pseudo = set(load_pseudo(cls))

    def ok(a):
        return a in pseudo and al.emit(a, cls) is not None

    cand = sorted(a for a in refcount if ok(a))
    rare = [a for a in cand if len(refcount[a]) <= rare_max]
    med = [a for a in cand if rare_max < len(refcount[a]) < freq_min]
    freq = [a for a in cand if len(refcount[a]) >= freq_min]
    rng.shuffle(med)
    rng.shuffle(freq)
    return sorted(set(rare) | set(med[:n_sample]) | set(freq[:n_sample]))


def train_records(refcount: dict, test: dict, cls: str) -> list:
    """mhcmatch training records with each held (epitope, allele) pair removed (per-pMHC)."""
    label = _LABEL[cls]
    return [{"epitope": p, "mhc_a": a, "mhc_class": label}
            for a, peps in refcount.items() for p in peps if p not in test.get(a, ())]


def holdout_split(refcount: dict, eval_alleles, cls: str, rng, frac: float = 0.3, cap: int = 40):
    """Per-pMHC holdout: hold out a capped fraction of each eval allele's ligands as positives.
    Returns ``(test={allele: set(held)}, train_records)``."""
    test = {}
    for a in sorted(eval_alleles):  # deterministic RNG-draw order (see select_eval_alleles)
        peps = list(refcount[a])
        rng.shuffle(peps)
        k = min(cap, max(1, int(frac * len(peps))))
        test[a] = set(peps[:k])
    return test, train_records(refcount, test, cls)


def loao_split(refcount: dict, held: str, cls: str):
    """Zero-shot leave-one-allele-out: ALL of ``held``'s ligands become positives and every one of
    them is removed from mhcmatch's training. Returns ``(test, train_records)``."""
    label = _LABEL[cls]
    test = {held: set(refcount[held])}
    train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
             for a, peps in refcount.items() if a != held for p in peps]
    return test, train


if __name__ == "__main__":
    import random

    rng = random.Random(0)
    rc = load_canonical(os.path.expanduser("~/hf/pmhc_data"), "mhc1", "human", tier="shortlist")
    assert all("*" not in a for a in rc), "class-I keys must be canonical (no star)"
    ev = select_eval_alleles(rc, "mhc1", rng, n_sample=10)
    rare = [a for a in ev if len(rc[a]) <= 30]
    assert ev and rare, (len(ev), len(rare))
    test, train = holdout_split(rc, ev, "mhc1", rng)
    held_pairs = {(a, p) for a, ps in test.items() for p in ps}
    train_pairs = {(r["mhc_a"], r["epitope"]) for r in train}
    assert held_pairs.isdisjoint(train_pairs), "held pairs leaked into training"
    a0 = rare[0]
    _, loao = loao_split(rc, a0, "mhc1")
    assert a0 not in {r["mhc_a"] for r in loao}, "LOAO allele still in training"
    print(f"splits.py self-check OK  {len(rc)} alleles, {len(ev)} eval ({len(rare)} rare); "
          f"holdout held {len(held_pairs)} pairs; LOAO removed {a0}")
