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
hydrophobic-run and dipeptide motifs, and per-role residue log-odds — whose two columns reproduce
:func:`mhcmatch.posbayes.llr` exactly, so that model is a strict special case. Linear head, because
a diagonal Gaussian cannot represent a summed log-odds; the EM Gaussian parameters ship alongside
for comparison. **Vectorised** — pass a list, not a loop. See :doc:`complementarity`.

.. automodule:: mhcmatch.complement
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
