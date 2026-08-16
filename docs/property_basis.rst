The amino-acid property basis
=============================

This page states two properties of the vendored tables in :mod:`mhcmatch.data.aa_tables`. They are
facts about the *data tables*, not about any peptide set, any label, or any benchmark — every number
below is reproducible from a fresh install with no download, by the snippet at the bottom of the
page. They are recorded here because both of them change what a physicochemical featurizer can
possibly express, and both are easy to get wrong.

1. The dominant eigenvector is a hydropathy axis.
2. The Kidera factors are already orthogonal, so PCA on them is degenerate.

What is vendored
----------------

``aa_tables`` carries 148 residue-to-value tables: 102 components across 17 descriptor families
(``DESCRIPTORS``), 45 hydrophobicity scales (``HYDROPHOBICITY``), and the Miyazawa-Jernigan
partition energy (``MJ_PARTITION``). 142 of the 148 are complete and non-constant over the standard
20 residues; the 6 that are not are the pH-variant Wimley-White scales, which are undefined for some
residues and are excluded from every statement here rather than imputed.

Provenance is in the module docstring: ``DESCRIPTORS`` and ``HYDROPHOBICITY`` are vendored from
peptides 0.5.0 (GPL-3.0-or-later), ``MJ_PARTITION`` from tcren 2.8.0, itself AAindex ``MIYS850101``.
The tables are generated, not hand-edited.

1. PC1 is hydropathy
--------------------

Column-standardize the 20 x 142 property matrix — **residues x scales**, not samples x features — and
take its principal components. The first carries **32.79 %** of the total variance and the second
**18.43 %**, so two axes hold just over half of everything the 142 scales measure. PC1's residue
order is

.. code-block:: text

   I F L W V M C Y A P G T H S Q N E K D R

which is a hydropathy ranking read off no hydropathy scale in particular. Against the 39 complete
named hydrophobicity scales, the Spearman rank correlation of PC1 has median :math:`|\rho| = 0.894`,
and 32 of the 39 reach :math:`|\rho| \ge 0.80`:

.. list-table:: Spearman correlation of PC1 with named hydrophobicity scales (n = 20 residues)
   :header-rows: 1

   * - scale
     - Spearman rho
     - p
   * - Fauchere
     - **+0.9774**
     - 1.3e-13
   * - Cowan7.5
     - **+0.9774**
     - 1.3e-13
   * - Eisenberg
     - +0.9700
     - 1.7e-12
   * - Roseman
     - +0.9643
     - 7.9e-12
   * - BlackMould
     - +0.9507
     - 1.4e-10
   * - Prabhakaran
     - -0.9383
     - 9.8e-10
   * - KyteDoolittle
     - +0.8943
     - 1.1e-07

(Bold marks the strongest correlation; Fauchere and Cowan7.5 tie exactly. The sign of a principal
component is arbitrary, so it is fixed here to make isoleucine positive; ``Prabhakaran`` is a
polarity scale and its negative sign is the same statement.)

PC2 is *not* a second reading of the same axis — against Kyte-Doolittle it sits at
:math:`\rho = +0.408` (p = 0.074), i.e. not distinguishable from unrelated at this sample size.

**What follows for the library.** Any summed or averaged physicochemical feature is, up to a
monotone reparametrization, a projection onto this one axis. So:

- Adding hydrophobicity scales to :func:`mhcmatch.immuno.features` adds *columns*, not *directions*;
  a disagreement between two such features is a disagreement about scaling, not about chemistry.
- A result that changes when the hydrophobicity scale is swapped is telling you something about the
  scale, not about the peptides.
- :mod:`mhcmatch.ipred` ships this basis rather than recomputing it: ``ipred.residue_scores()``
  returns the frozen per-residue PC1/PC2 loadings, and its ``pc1`` feature is the sum of PC1 along
  the sequence.

2. PCA on the Kidera factors is a no-op
---------------------------------------

The 10 Kidera factors are the orthogonal factor scores of a 1985 factor analysis of 188 physical
properties. Orthogonality is not an approximation that decayed in transcription — it survives in the
vendored table exactly. Over the 20 residues the largest absolute off-diagonal correlation between
two Kidera factors is **0.0026**, and the correlation matrix's eigenvalues are

.. code-block:: text

   1.0043 1.0039 1.0024 1.0014 1.0004 0.9998 0.9988 0.9976 0.9965 0.9949

— a ratio of 1.0094 from largest to smallest, and a participation ratio of 10.00 out of 10
components. A covariance matrix whose eigenvalues are all equal is a sphere: its eigenbasis is
defined only up to an arbitrary rotation, and PCA on it returns a rotation and reduces nothing.

.. list-table:: Isotropy of three property bases over the 20 residues, uniform measure
   :header-rows: 1

   * - table
     - components
     - rank
     - max abs. off-diagonal correlation
     - largest / smallest non-zero eigenvalue
     - participation ratio
   * - ``DESCRIPTORS["KIDERA"]``
     - 10
     - 10
     - 0.0026
     - 1.01
     - 10.00
   * - ``DESCRIPTORS["VHSE"]``
     - 8
     - 8
     - 0.8201
     - 50.43
     - 4.79
   * - all 142 complete scales
     - 142
     - 19
     - 1.0000
     - 101.23
     - 5.99

(Nothing is bolded: this is a description of three bases, not a competition between them. A
participation ratio equal to the rank means perfectly isotropic — Kidera is at 10.00 of 10 — and the
further below the rank it falls, the more the variance concentrates in a few directions. The union
of all 142 scales has rank 19, so 123 of those columns are exactly linearly dependent on the others,
and its largest off-diagonal correlation of 1.0000 says at least one pair is the same scale twice
under different names.)

**The rotation destroys the hydropathy axis rather than finding it.** Because the sphere has no
preferred direction, the PC1 that a solver happens to return from the Kidera table alone correlates
with Kyte-Doolittle at only :math:`\rho = -0.256` (p = 0.28). Hydropathy is not spread across the
Kidera factors — it sits almost entirely in KF4:

.. list-table:: Each Kidera factor against Kyte-Doolittle (n = 20 residues)
   :header-rows: 1

   * - factor
     - Spearman rho
     - p
   * - KF1
     - -0.3042
     - 0.19
   * - KF2
     - -0.4053
     - 0.076
   * - KF3
     - +0.2348
     - 0.32
   * - KF4
     - **-0.7761**
     - 5.8e-05
   * - KF5
     - +0.0370
     - 0.88
   * - KF6
     - -0.2635
     - 0.26
   * - KF7
     - -0.0355
     - 0.88
   * - KF8
     - +0.0298
     - 0.90
   * - KF9
     - +0.0015
     - 0.99
   * - KF10
     - -0.0468
     - 0.84

So: use the raw factors, or run PCA over the multi-family union where the collinearity is real.
Running it over Kidera alone spends compute to scramble an interpretable basis.

Scope: what the degeneracy is *not*
-----------------------------------

The no-op holds for the 20 x 10 alphabet table under the uniform measure over the 20 residue types.
It is not a general licence to skip PCA anywhere Kidera factors appear:

- **Reweight the alphabet and it breaks.** Under residue frequencies rather than one weight per
  residue type — for example :data:`mhcmatch.diffusion.PROTEOME_AA_FREQ` — the same table gives an
  eigenvalue ratio of 3.72 and a participation ratio of 8.92 of 10, not 10.00 of 10. Isotropy is a
  property of the table *and* the measure.
- **Sum along a peptide and it breaks.** A ``peptides x 10 summed-Kidera`` matrix inherits the
  composition statistics of the peptide set and is strongly anisotropic; PCA on that matrix does
  real work. Anyone reading "PCA on Kidera is degenerate" as covering that case will be wrong.

Reproduce
---------

No data download, no benchmark repo, numpy only:

.. code-block:: python

   import numpy as np
   from mhcmatch.data import aa_tables as t

   AA = list(t.AA20)

   def spectrum(tables):
       """Eigenvalues of the residue-by-scale correlation matrix, uniform over residues."""
       A = np.array([[tab[a] for tab in tables] for a in AA])
       Z = (A - A.mean(0)) / A.std(0, ddof=0)
       ev = np.linalg.eigvalsh(Z.T @ Z / len(AA))[::-1]
       return ev[ev > 1e-9]

   kidera = [t.DESCRIPTORS["KIDERA"][f"KF{i}"] for i in range(1, 11)]
   ev = spectrum(kidera)
   ev[0] / ev[-1]                          # 1.0094 -- a sphere
   ev.sum() ** 2 / (ev ** 2).sum()         # 10.00 of 10 components

   scales = [tab for fam in t.DESCRIPTORS.values() for tab in fam.values()]
   scales += list(t.HYDROPHOBICITY.values()) + [t.MJ_PARTITION]
   scales = [s for s in scales if all(a in s for a in AA) and len(set(s.values())) > 1]
   len(scales)                             # 142

   A = np.array([[s[a] for s in scales] for a in AA])
   Z = (A - A.mean(0)) / A.std(0, ddof=0)
   U, S, Vt = np.linalg.svd(Z / np.sqrt(len(AA)), full_matrices=False)
   (S ** 2 / (S ** 2).sum())[:2]           # 0.3279, 0.1843

   pc1 = Z @ Vt[0]
   "".join(np.array(AA)[np.argsort(-pc1)])  # IFLWVMCYAPGTHSQNEKDR (up to overall sign)

The two claims are asserted as regression tests in ``tests/test_mhcmatch.py``
(``test_property_matrix_pc1_is_a_hydropathy_axis``, ``test_kidera_table_is_orthogonal``,
``test_kidera_degeneracy_is_specific_to_the_uniform_alphabet_measure``), so a table regenerated with
different provenance cannot silently change them.
