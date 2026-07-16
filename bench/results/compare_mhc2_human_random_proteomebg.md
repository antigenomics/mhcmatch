# mhcmatch vs NetMHCIIpan-4.3i (holdout, random decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 19 | **0.884** | 0.881 | +0.003 | [-0.075, 0.083] | 0.861 |
| rare | AUPRC | 19 | **0.652** | 0.610 | +0.042 | [-0.150, 0.225] | 0.625 |
| rare | PPV@P | 19 | **0.541** | 0.518 | +0.023 | [-0.233, 0.279] | 0.881 |
| medium | AUROC | 8 | 0.826 | **0.894** | -0.068 | [-0.127, -0.018] | 0.000 |
| medium | AUPRC | 8 | 0.460 | **0.574** | -0.114 | [-0.206, -0.032] | 0.000 |
| medium | PPV@P | 8 | 0.421 | **0.556** | -0.135 | [-0.231, -0.050] | 0.000 |
| frequent | AUROC | 20 | 0.884 | **0.966** | -0.082 | [-0.124, -0.044] | 0.000 |
| frequent | AUPRC | 20 | 0.524 | **0.775** | -0.250 | [-0.358, -0.154] | 0.000 |
| frequent | PPV@P | 20 | 0.496 | **0.733** | -0.236 | [-0.316, -0.160] | 0.000 |

## Re-baseline: `register="marginal"` is now the default (v0.6)

The table above is `--register marginal`; the previous default was `max` (max over 9-mer frames).
Identical examples, NetMHC scores, and seed — only `AnchorModel.score` changed
(`run_compare.py --register {max,marginal}`). See `register_em_mhc2.md` for the mechanism.

| stratum | metric | mhcmatch `max` (old default) | mhcmatch `marginal` (new) | Δ |
|---|---|---|---|---|
| rare | AUROC | 0.866 | **0.884** | +0.018 |
| rare | AUPRC | 0.555 | **0.652** | +0.097 |
| rare | PPV@P | 0.376 | **0.541** | +0.165 |
| medium | AUROC | 0.810 | **0.826** | +0.016 |
| medium | AUPRC | 0.446 | **0.460** | +0.014 |
| medium | PPV@P | 0.383 | **0.421** | +0.038 |
| frequent | AUROC | 0.874 | **0.884** | +0.010 |
| frequent | AUPRC | 0.467 | **0.524** | +0.057 |
| frequent | PPV@P | 0.451 | **0.496** | +0.045 |

**Every stratum × metric improves; none regresses.** The gain is largest on **rare** alleles
(AUPRC +0.097, PPV@P +0.165), which flip from losing all three metrics to winning all three — not
significant at n=19 (p=0.86 / 0.63 / 0.88), but the direction is consistent across both decoy modes.
That ordering is mechanistic: where the motif is too thin to pin the register on its own, the offset
prior carries proportionally more of the decision.

**The screening gap remains the real one.** Frequent AUPRC closes -0.308 → **-0.250** (19% of it) and
NetMHCIIpan still leads decisively. This is the task where a presentation background was already the
headline fix (`ROADMAP.md` §6b) and it is still not enough — consistent with the gap being capacity
or memorisation rather than register handling. Note NetMHCIIpan is run **without** `-context`
(`bench/compare/netmhc.py:95`, `-inptype 1` peptide list), so this gap is not a flank-model gap.
