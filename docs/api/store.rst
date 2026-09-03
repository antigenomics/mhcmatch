Store, search and reference data
--------------------------------

The entry points, and the deposits everything else reads from.

mhcmatch.store module
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mhcmatch.store
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.search module
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mhcmatch.search
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.proteome module
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mhcmatch.proteome
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.pseudoseq module
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mhcmatch.pseudoseq
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.expression module
~~~~~~~~~~~~~~~~~~~~~~~~~~

Reference expression by normal tissue (GTEx) and by tumour type (TCGA), fetched from the public
``isalgo/pmhc_data`` dataset. The two are never merged — different measurements, different units.

.. automodule:: mhcmatch.expression
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.known module
~~~~~~~~~~~~~~~~~~~~~

Built-in known-epitope reference sets for exact-match lookup, assembled from the public deposits:
confirmed tumour neoantigens, peptides the screens tested and found **negative**, IEDB-immunogenic
epitopes, the thymic self-immunopeptidome and the viral ligandome. An exact match is stronger
evidence than any model output, so :mod:`mhcmatch.rank` reports it as a flag and never folds it
into the score.

.. automodule:: mhcmatch.known
   :members:
   :undoc-members:
   :show-inheritance:
