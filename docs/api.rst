API reference
=============

mhcmatch.store module
---------------------

.. automodule:: mhcmatch.store
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.search module
----------------------

.. automodule:: mhcmatch.search
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.proteome module
------------------------

.. automodule:: mhcmatch.proteome
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.pseudoseq module
-------------------------

.. automodule:: mhcmatch.pseudoseq
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.diffusion module
-------------------------

.. automodule:: mhcmatch.diffusion
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.calibrate module
-------------------------

.. automodule:: mhcmatch.calibrate
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.affinity module
------------------------

.. automodule:: mhcmatch.affinity
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.structure module
-------------------------

.. automodule:: mhcmatch.structure
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.ligand module
----------------------

.. automodule:: mhcmatch.ligand
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.logo module
--------------------

.. automodule:: mhcmatch.logo
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.predict module
-----------------------

.. automodule:: mhcmatch.predict
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.mimics module
----------------------

Cross-reactivity by **category**, never summed: a hit in the thymic immunopeptidome argues tolerance
and autoimmune risk, a hit in the host proteome argues peripheral cross-reactivity without implying
presentation, and a viral or bacterial hit argues a pre-existing repertoire may cross-react.
:data:`~mhcmatch.mimics.KINDS` states which is which; :func:`~mhcmatch.mimics.neighbours` runs the
whole query set through one threaded C++ index per (category, length).

.. automodule:: mhcmatch.mimics
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.immuno module
----------------------

.. automodule:: mhcmatch.immuno
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.ipred module
---------------------

.. automodule:: mhcmatch.ipred
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.posbayes module
------------------------

Position-role naive Bayes over amino-acid identity: anchor and TCR-facing residues get separate
conditional distributions, because for several amino acids their contributions carry opposite signs
and pooling averages that away. Emits a **log-likelihood ratio**, so the caller supplies the prior.

Grouped 5-fold CV **0.712** human / **0.758** mouse (against ``ipred`` *in-sample* at 0.607 / 0.668);
size-matched transfer **0.731** human→mouse, **0.692** mouse→human. Cysteine is masked — see the
module warning.

.. automodule:: mhcmatch.posbayes
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.complement module
--------------------------

The recognition axis as one score: ``ipred``'s physicochemistry and length, the same components
split MHC-facing vs TCR-facing, MJ1996 / repertoire-marginalised TCRen contact potentials,
hydrophobic-run and dipeptide motifs, and per-role residue log-odds — pooled, **per length bin
(8/9/10/11+)** and per relative third of the TCR face, whose pooled ``aa_anchor``/``aa_tcr`` pair
reproduces :func:`mhcmatch.posbayes.llr` exactly, so that model is a strict special case. Linear
head, because
a diagonal Gaussian cannot represent a summed log-odds; the EM Gaussian parameters ship alongside
for comparison. **Vectorised** — pass a list, not a loop. See :doc:`complementarity`.

.. automodule:: mhcmatch.complement
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.mimicry module
-----------------------

Mimicry as a signed, per-component immune-response risk. Three references — ``viral`` (priming),
``self`` (tolerance, and the autoimmunity read-out) and ``thymus`` (negative selection) — each split
into an **anchor** and a **TCR-facing** channel, because a whole-peptide distance averages two
different measurements. Scores are **log-odds**; :func:`mhcmatch.mimicry.probability` is a separate
step that requires a *named* corpus, because the screens behind any calibration run from 0.048 % to
46.8 % positive and an unqualified probability is mostly a statement about which intercept was used.

The tested-neoantigen database is exposed as :func:`mhcmatch.mimicry.annotate` — **prior evidence,
never a fitted term**. Every labelled screen we hold sits inside that database, so a coefficient on
it would be memorisation; held out honestly, fuzzy matching at two substitutions still recovers
0.08–0.34 of a fresh screen's positives against 0.00–0.26 for exact lookup, which is what makes it
worth reporting. Command line: ``mhcmatch neoag``.

.. automodule:: mhcmatch.mimicry
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.known module
---------------------

Built-in known-epitope reference sets for exact-match lookup, assembled from the public deposits:
confirmed tumour neoantigens, peptides the screens tested and found **negative**, IEDB-immunogenic
epitopes, the thymic self-immunopeptidome and the viral ligandome. An exact match is stronger
evidence than any model output, so :mod:`mhcmatch.rank` reports it as a flag and never folds it
into the score.

.. automodule:: mhcmatch.known
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.expression module
--------------------------

Reference expression by normal tissue (GTEx) and by tumour type (TCGA), fetched from the public
``isalgo/pmhc_data`` dataset. The two are never merged — different measurements, different units.

.. automodule:: mhcmatch.expression
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.rank module
--------------------

Neoantigen candidate ranking, from a mutation-spanning window FASTA or an already-scored table.
Combines presentation and recognition through a **gate** (a product of sigmoids) rather than a sum.

.. automodule:: mhcmatch.rank
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.precursor module
-------------------------

Optional extra: ``pip install 'mhcmatch[precursor]'`` (needs ``vdjtools``).

.. automodule:: mhcmatch.precursor
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.data.aa_tables module
------------------------------

The vendored amino-acid property tables. Their values are not reproduced here — they are plain
``dict[str, float]`` and large; :doc:`property_basis` states what their principal components are.

.. automodule:: mhcmatch.data.aa_tables

mhcmatch.data.contact_profile module
------------------------------------

.. automodule:: mhcmatch.data.contact_profile
