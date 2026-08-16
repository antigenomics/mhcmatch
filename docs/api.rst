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
