# mhcmatch vs NetMHCpan — leak-free affinity on TESLA measured IC50

Both predictors scored against TESLA-608's **MEASURED_BINDING_AFFINITY** (competition binding), the held-out ground truth that avoids the IEDB train/test overlap inflating NetMHCpan in `bench/affinity/eval.py`. **496 candidates** with a measured value (363 binders ≤500 nM). NetMHCpan prediction is the dataset's embedded column (zero rerun). **Bold = better.**

| metric | mhcmatch | NetMHCpan |
|---|--:|--:|
| Spearman(pred, measured) | 0.628 | **0.650** |
| AUROC @500 nM | 0.799 | **0.813** |

Per-allele (macro over 10 alleles with ≥8 measured points): median Spearman mhcmatch 0.712 vs NetMHCpan 0.684; median AUROC 0.844 vs 0.827.

> Read with the training caveat: NetMHCpan may have seen some TESLA-adjacent IEDB measurements, mhcmatch's Potts was fit on the complementary IEDB split; TESLA is closer to out-of-sample for both than the in-corpus `eval.py` holdout.
