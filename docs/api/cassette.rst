Cassette design
---------------

Selecting units, assembling a construct, and what a portfolio is worth.

mhcmatch.cassette module
~~~~~~~~~~~~~~~~~~~~~~~~

Choosing the units of a cassette, and scoring one that already exists. The narrative version, with
the derivation and the measured numbers, is :doc:`cassette`.

.. automodule:: mhcmatch.cassette
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.vector module
~~~~~~~~~~~~~~~~~~~~~~

Assembling the cassette, once selection has chosen the candidates: **what to withdraw on safety
grounds, how many units each allotype should carry, in what order, and joined by what.**

:func:`~mhcmatch.vector.screen` runs first, and it excludes rather than down-ranks: a register that
is itself an essential-tissue self peptide, or a target gene transcribed where it was assumed silent,
is withdrawn, because the second-best cassette is cheap and myocarditis is not. The two fatal
precedents behind that rule, and the measurement that chose
:func:`~mhcmatch.vector.self_origin_risk` over a mimicry-similarity screen, are in the function's own
documentation. Competition for a response is then local to the
antigen-presenting cell and strongest *within* an allotype, so expected yield is a sum of
independently saturating per-allotype terms rather than one global budget --
:func:`~mhcmatch.vector.select` grows each allotype while the next candidate beats that allotype's
own expected yield per slot, and diversification follows from the arithmetic instead of a quota.
:func:`~mhcmatch.vector.order` scores every register spanning every junction against the recipient's
own allotypes and picks the spacer and ordering that minimise predicted junctional binding, trying
**no spacer first**.

The two ends of that pipeline are joins to the rest of the library.
:func:`~mhcmatch.vector.units_from_context` closes the front: :func:`~mhcmatch.rank.rank_fasta`
emits *minimal epitopes* while a unit is the long window around the mutation, and where that
mutation sits is in the FASTA header rather than in the ranking, so the two are combined rather than
either being guessed at -- and rows are grouped by **variant**, since twenty registers of one
mutation are twenty rows in a ranking and one thing to put in a cassette.
:func:`~mhcmatch.vector.back_translate` closes the back, turning
:attr:`~mhcmatch.vector.Cassette.sequence` into a coding sequence. It is not a codon optimiser: it
fixes the two failure modes specific to a *concatemer* -- the m1-pseudouridine +1-frameshift motif
(:func:`~mhcmatch.vector.slippery_sites`, whose seams the designer chooses) and synthesis-hostile
homopolymers (which spacers like ``AAA`` manufacture directly) -- and leaves GC content, secondary
structure and CpG to a manufacturer's own tooling. :func:`~mhcmatch.vector.translate` exists so
"synonymous" stays checkable rather than asserted.

.. automodule:: mhcmatch.vector
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.portfolio module
~~~~~~~~~~~~~~~~~~~~~~~~~

The composition layer above :func:`~mhcmatch.vector.select`: objective-space geometry (Pareto front,
crowding, hull membership, Chebyshev scalarization) and the block response model that says what a
proposed cassette is worth. Narrative and worked examples: :doc:`../portfolio`.

.. automodule:: mhcmatch.portfolio
   :members:
   :undoc-members:
   :show-inheritance:
