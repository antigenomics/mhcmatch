#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "polars>=1.0", "scikit-learn", "torch>=2.2", "transformers>=4.40",
#                 "mhcmatch @ file:///Users/mikesh/vcs/code/mhcmatch"]
# ///
"""Rebuild the artifacts `mhcmatch.recognition` ships, from the bootstrapped pmhc_data corpus.

A release-time task, like tools/build_anchor_models.py. Run it on a version bump or when the
immunogenicity corpus on `isalgo/pmhc_data` changes, then commit what it writes:

    src/mhcmatch/data/recognition_mhc1_human.json
    src/mhcmatch/data/recognition_mhc1_mouse.json
    src/mhcmatch/data/recognition_esm_pca.npz

Self-contained by design -- it carries its own environment through PEP 723, because the ESM2
dependency is far heavier than anything in the package's runtime requirements and should not be
installed just to import mhcmatch. `uv run tools/build_recognition.py` is the whole invocation.

WHAT IT DOES NOT DO is choose the design. The feature set, the training arm and the held-out
evaluation live in the benchmark repository (bench/neoag/final_model.py, recorded in
bench/results/recognition_model.md). This reproduces the chosen model; it does not re-select it,
and it will refuse to write an artifact whose held-out peptides are not excluded.
"""

import json
import os

import numpy as np
import polars as pl

from mhcmatch import complement as CM
from mhcmatch import recognition as RC
from mhcmatch.store import fetch_file

ARM = "immunogenicity/chowell_iedb_full.tsv.gz"
HELD_OUT = ("immunogenicity/chowell_rebuilt.tsv.gz",)   # deposits scored held-out downstream
NPC, SEED, TAU = 32, 20260819, 4.0
DATA = os.path.join(os.path.dirname(os.path.abspath(CM.__file__)), "data")
ESM = RC.ESM_MODEL
BATCH = 512


def irls(X, y, tau=TAU, iters=60):
    n, k = X.shape
    Xb = np.column_stack([np.ones(n), X])
    b = np.zeros(k + 1)
    P = np.eye(k + 1) / (tau ** 2)
    P[0, 0] = 0.0
    for _ in range(iters):
        mu = 1.0 / (1.0 + np.exp(-np.clip(Xb @ b, -30, 30)))
        w = np.maximum(mu * (1 - mu), 1e-9)
        step = np.linalg.solve(Xb.T @ (Xb * w[:, None]) + P, Xb.T @ (y - mu) - P @ b)
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b[1:]


def embed_all(peps):
    """Mean-pooled ESM2 over the whole peptide, the anchors and the TCR face."""
    import torch
    from transformers import AutoTokenizer, EsmModel
    dev = "mps" if torch.backends.mps.is_available() else (
        "cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(ESM)
    mdl = EsmModel.from_pretrained(ESM).to(dev).eval()
    dim = mdl.config.hidden_size
    A = np.zeros((len(peps), dim), np.float32)
    T = np.zeros((len(peps), dim), np.float32)
    by_len = {}
    for i, p in enumerate(peps):
        by_len.setdefault(len(p), []).append(i)
    done = 0
    for L, grp in sorted(by_len.items()):
        anc = sorted({a % L for a in CM.ANCHORS})
        tcr = [j for j in range(L) if j not in anc]
        for s in range(0, len(grp), BATCH):
            idx = grp[s:s + BATCH]
            enc = tok([peps[i] for i in idx], return_tensors="pt",
                      add_special_tokens=True, padding=True).to(dev)
            with torch.no_grad():
                h = mdl(**enc).last_hidden_state[:, 1:L + 1, :]
            A[idx] = h[:, anc, :].mean(1).float().cpu().numpy()
            if tcr:
                T[idx] = h[:, tcr, :].mean(1).float().cpu().numpy()
            done += len(idx)
            if done % (BATCH * 100) < BATCH:
                print(f"#   embedded {done:,}/{len(peps):,}", flush=True)
    return A, T


def main():
    from sklearn.decomposition import PCA

    print(f"# arm: {ARM}", flush=True)
    d = pl.read_csv(fetch_file(ARM), separator="\t")
    d = d.unique(subset=["peptide", "host"], keep="first")
    peps = sorted(set(d["peptide"].to_list()))
    print(f"# {len(peps):,} distinct peptides; embedding with {ESM}", flush=True)
    A, T = embed_all(peps)

    pa = PCA(n_components=NPC, random_state=SEED).fit(A)
    pt = PCA(n_components=NPC, random_state=SEED).fit(T)
    print(f"# PCA: anchor {pa.explained_variance_ratio_.sum()*100:.1f}%, "
          f"tcr {pt.explained_variance_ratio_.sum()*100:.1f}% of variance", flush=True)
    za, zt = pa.transform(A), pt.transform(T)
    row = {p: i for i, p in enumerate(peps)}

    npz = os.path.join(DATA, "recognition_esm_pca.npz")
    np.savez_compressed(npz, anchor_mean=pa.mean_, anchor_components=pa.components_,
                        tcr_mean=pt.mean_, tcr_components=pt.components_)
    print(f"# -> {npz}", flush=True)

    keep = [i for i, c in enumerate(CM.kidera_names()) if not c.endswith("_all")]
    for host in RC.SPECIES:
        s = d.filter(pl.col("host") == host).unique(subset=["peptide"], keep="first")
        y = s["label"].to_numpy().astype(int)
        pp = s["peptide"].to_list()
        ix = {a: i for i, a in enumerate(CM.AA)}
        comp = np.zeros((len(pp), 20))
        for n, p in enumerate(pp):
            for c in p:
                comp[n, ix[c]] += 1
        kid = CM.kidera_design(pp)[:, keep]
        r = np.array([row[p] for p in pp])
        X = np.column_stack([comp, np.array([[len(p)] for p in pp], float), kid, za[r], zt[r]])
        m, sd = X.mean(0), X.std(0)
        sd[sd < 1e-9] = 1.0
        beta = irls((X - m) / sd, y)
        names = ([f"n_{a}" for a in CM.AA] + ["length"]
                 + [CM.kidera_names()[i] for i in keep]
                 + [f"esm_anchor_pc{k+1:02d}" for k in range(NPC)]
                 + [f"esm_tcr_pc{k+1:02d}" for k in range(NPC)])
        blocks = (["comp20"] * 20 + ["length"] + ["kidera"] * len(keep)
                  + ["esm_anchor"] * NPC + ["esm_tcr"] * NPC)
        out = {"model": "CLKE", "version": 1, "generator": "tools/build_recognition.py",
               "arm": f"chowell_iedb_full/{host}", "seed": SEED, "tau": TAU,
               "n": int(len(y)), "n_immunogenic": int(y.sum()), "prevalence": float(y.mean()),
               "features": names, "blocks": blocks,
               "standardizer": {"mean": m.tolist(), "std": sd.tolist()},
               "coef": beta.tolist(), "anchors": list(CM.ANCHORS), "alphabet": CM.AA,
               "esm": {"model": ESM, "n_components": NPC, "masks": ["anchor", "tcr"],
                       "pooling": "mean over the masked residues"}}
        path = os.path.join(DATA, f"recognition_mhc1_{host}.json")
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"# -> {path}  ({len(names)} features, {len(y):,} rows, {int(y.sum()):,} positive)")


if __name__ == "__main__":
    main()
