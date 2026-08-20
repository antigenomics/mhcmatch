Physicochemical immunogenicity features
=======================================

**What this shows.** The path from a fresh install to a physicochemical feature vector for your own
peptides, and what the dominant axis of that vector is.

**What you should conclude.** :func:`mhcmatch.immuno.features` is not a bag of descriptors with two
free choices left implicit. *Which* positions count as TCR-facing and *how* residues are aggregated
are both arguments you pass, because neither is settled; and the scales themselves are so collinear
that a single hydropathy axis carries a third of their total variance
(:doc:`property_basis`).

Nothing on this page touches the reference panel. The scales and the contact profile are vendored in
the package, so every snippet runs offline in under a second and needs no HuggingFace download.

Install
-------

.. code-block:: bash

   pip install mhcmatch

   python -m mhcmatch.immuno      # self-check: prints "ok - 141 features, 20 scales, ..."

One peptide, 141 numbers
------------------------

.. code-block:: python

   from mhcmatch import immuno

   f = immuno.features("GILGFVFTL")        # influenza A M1 58-66, HLA-A*02:01
   len(f)                                  # 141
   len(immuno.scales())                    # 20  ->  1 length + 20 scales x 7 statistics
   immuno.feature_names()[:8]
   # ['length', 'KF1_sum', 'KF1_mean', 'KF1_min', 'KF1_max',
   #  'KF1_run_max', 'KF1_run_n', 'KF1_run_frac']

The seven statistics per scale are ``sum``, ``mean``, ``min``, ``max``, ``run_max``, ``run_n``,
``run_frac``. The first four are the established descriptors. The three ``run_*`` ones exist because
a *contiguous* hydrophobic stretch is a different object from the same residues scattered along the
peptide, and no sum can express that — a permutation of a peptide has, by construction, identical
``sum``/``mean``/``min``/``max``, and different ``run_max``/``run_n``.

The default scale set is :data:`mhcmatch.immuno.DEFAULT_SCALES` — the 10 Kidera factors, the 8 VHSE
components, the Miyazawa-Jernigan partition energy, and Kyte-Doolittle. Any table in
:mod:`mhcmatch.data.aa_tables` can be substituted::

   immuno.features("GILGFVFTL", scale_names=("KyteDoolittle", "Eisenberg", "MJ"))

Which positions count as TCR-facing
-----------------------------------

Three incompatible class-I anchor definitions coexist in this toolchain, so all three are kept
selectable in :data:`mhcmatch.immuno.ANCHOR_SCHEMES` rather than collapsed into a constant. A fourth
option, ``"contact"``, needs no anchor call at all: it weights each position by observed
TCR-peptide contact frequency over 370 crystal structures.

.. code-block:: python

   immuno.position_weights("GILGFVFTL", "mhc1", "p2_pomega")
   # [1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]        masks P2 + POmega

   cp = immuno.contact_profile("mhc1")
   [round(w, 2) for w in cp(9)]
   # [0.0, 0.0, 0.0, 1.26, 1.14, 1.03, 0.7, 0.86, 0.0]     masks P1, P2, P3, POmega

   immuno.features("GILGFVFTL", scheme="contact", contact_profile=cp)

The scheme changes the numbers, and it is meant to — it is an ablation axis with a reported number,
not a default to accept silently:

.. list-table:: Kyte-Doolittle statistics for ``GILGFVFTL`` under each scheme
   :header-rows: 1

   * - scheme
     - masked positions
     - ``KyteDoolittle_sum``
     - ``KyteDoolittle_mean``
     - ``KyteDoolittle_run_max``
     - ``KyteDoolittle_run_n``
   * - ``full``
     - 0
     - 20.400
     - 2.267
     - 9
     - 1
   * - ``p2_pomega``
     - 2
     - 12.100
     - 1.729
     - 6
     - 2
   * - ``pockets``
     - 5
     - 9.400
     - 2.350
     - 4
     - 1
   * - ``contact``
     - 4
     - 8.394
     - 1.679
     - 5
     - 1

No cell is bolded because no scheme is the winner: which one to use is the question, and answering
it is a benchmark, not a library default. Note ``run_n`` under ``p2_pomega``: masking P2 *breaks* the
run rather than bridging it, so one stretch of 9 becomes two. A buried anchor between two exposed
hydrophobics does not make them contiguous from the TCR's point of view.

Class II ignores ``scheme`` entirely and always masks the register-anchored core P1/P4/P6/P9,
because *that* definition is agreed across the toolchain. Pass ``register_start=`` from
:meth:`mhcmatch.diffusion.AnchorModel.best_register` so the annotated frame matches the scored one.

A feature matrix
----------------

:func:`mhcmatch.immuno.feature_names` gives the column order without a dict round-trip, so a matrix
is one comprehension:

.. code-block:: python

   peptides = ["GILGFVFTL", "NLVPMVATV", "GLCTLVAML", "TPRVTGGGAM"]
   cp = immuno.contact_profile("mhc1")

   cols = immuno.feature_names()
   X = [[immuno.features(p, scheme="contact", contact_profile=cp)[c] for c in cols]
        for p in peptides]                       # 4 x 141

Two conventions worth knowing before you fit anything on ``X``:

- **``length`` is a feature by decision, not an oversight.** The length distribution of an allele's
  ligand set is part of what defines it, so it is signal here rather than a nuisance to regress out.
- **Non-standard residues are dropped, not zeroed.** ``0`` is a real value on a centred scale such as
  Kidera, so scoring ``X`` as zero would be a silent bias.

What the dominant axis means
----------------------------

The 142 complete scales in :mod:`mhcmatch.data.aa_tables` are massively collinear. The first
principal component of the column-standardized 20 x 142 *property* matrix (residues x scales) carries
**32.79 %** of their total variance, and it is a hydropathy axis — its residue order is

.. code-block:: text

   I F L W V M C Y A P G T H S Q N E K D R

so a *summed* physicochemical feature is, up to a monotone reparametrization, a projection of the
peptide onto that one axis. This has a practical consequence: adding more hydrophobicity scales to
``scale_names`` adds columns but almost no directions. :doc:`property_basis` states the measurement,
including why running PCA on the Kidera factors specifically is a no-op.

The axis itself is shipped, so you do not have to recompute it — the frozen per-residue loadings
are :data:`mhcmatch.data.aa_tables.PROPERTY_PC1` and
:data:`~mhcmatch.data.aa_tables.PROPERTY_PC2`:

.. code-block:: python

   from mhcmatch.data import aa_tables

   sorted(aa_tables.PROPERTY_PC1, key=lambda a: -aa_tables.PROPERTY_PC1[a])
   # ['I', 'F', 'L', 'W', 'V', 'M', 'C', 'Y', 'A', 'P', 'G', 'T', 'H', 'S', 'Q', 'N', 'E', 'K', 'D', 'R']

The fitted model built on that axis is :func:`mhcmatch.complement.score`, which reads PC1/PC2 and
length as its ``phys`` block and adds five more (:doc:`complementarity`). Its predecessor
:mod:`!mhcmatch.ipred` — two components plus length, thirteen parameters, a calibrated
``P(immunogenic)`` — **shipped through 0.21.0 and was removed in 0.22.0**; what it was, how it
scored and why it went is recorded in full at :ref:`ipred-legacy`.

Use :func:`mhcmatch.complement.score` when you want the shipped answer and
:func:`mhcmatch.immuno.features` when you want to fit your own. Neither is an alternative to the
presentation heads: presentation asks whether a peptide reaches the surface, this asks which of the
peptides that do are recognised.

Further worked examples
-----------------------

``notebooks/02_immunogenicity_features.py`` (a marimo notebook, ``pip install 'mhcmatch[notebooks]'``)
runs the same material interactively and adds two demonstrations that need live output: three
peptides with identical amino-acid composition separated only by the run statistics, and the contact
profile recovering the class-I anchor set unsupervised.
