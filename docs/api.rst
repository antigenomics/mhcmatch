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

mhcmatch.vector module
----------------------

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
-------------------------

The composition layer above :func:`~mhcmatch.vector.select`: objective-space geometry (Pareto front,
crowding, hull membership, Chebyshev scalarization) and the block response model that says what a
proposed cassette is worth. Narrative and worked examples: :doc:`portfolio`.

.. automodule:: mhcmatch.portfolio
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.immuno module
----------------------

.. automodule:: mhcmatch.immuno
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.posbayes module
------------------------

Position-role naive Bayes over amino-acid identity: anchor and TCR-facing residues get separate
conditional distributions, because for several amino acids their contributions carry opposite signs
and pooling averages that away. Emits a **log-likelihood ratio**, so the caller supplies the prior.

Grouped 5-fold CV **0.712** human / **0.758** mouse (against the retired ``ipred`` *in-sample* at
0.607 / 0.668 --- :ref:`ipred-legacy`);
size-matched transfer **0.731** human→mouse, **0.692** mouse→human. Cysteine is masked — see the
module warning.

.. automodule:: mhcmatch.posbayes
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.complement module
--------------------------

The recognition axis as one score: the retired ``ipred``'s physicochemistry and length, the same components
split MHC-facing vs TCR-facing, MJ1996 / repertoire-marginalised TCRen contact potentials,
hydrophobic-run and dipeptide motifs, and per-role residue log-odds — pooled, **per length bin
(8/9/10/11+ at class I, quartiles at 14/16/19 at class II)** and per position zone (relative thirds
of the TCR face at class I, the ``nflank``/``core``/``cflank`` register zones at class II, selected
by ``cls=``), whose pooled ``aa_anchor``/``aa_tcr`` pair
reproduces :func:`mhcmatch.posbayes.llr` exactly, so that model is a strict special case. Linear
head, because
a diagonal Gaussian cannot represent a summed log-odds; the EM Gaussian parameters ship alongside
for comparison. **Vectorised** — pass a list, not a loop. See :doc:`complementarity`.

.. automodule:: mhcmatch.complement
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.recognition module
---------------------------

The head dispatcher over the recognition axis: ``complement`` (the six-block model above, and the
default), plus ``posbayes``, ``physchem_glm`` and ``esm64_glm``, each fitted alone so their BIC is
comparable. See :doc:`complementarity`.

.. automodule:: mhcmatch.recognition
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.luksza module
----------------------

The Łuksza recognition term :math:`R = Z/(1+Z)` -- a soft partition function over near-matches,
replacing a hard distance cut, so **how many** near-matches a candidate has and **how near** they
are both enter. ``viral_R`` was a term of the retired ``BOECRT`` aggregate — before 0.17.0 the sum
lived only in the benchmark repo, which made :func:`mhcmatch.rank.aggregate_score` callable with a
feature no installed user could supply. ``EPIC`` retired the term in 0.21.0, so the shipped
aggregate no longer scores with it; the quantity is still the published recognition term and is
still computed here. ``k`` and ``a0`` are read from the shipped artifact when one carries them, so
a refit needs no code change.

.. automodule:: mhcmatch.luksza
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
The default score is the fitted aggregate vendored at ``mhcmatch/data/aggregate_mhc1.json``
(``--score aggregate``); the noisy-AND **gate** — a product of sigmoids, the default before
0.19.0 — is still reachable as ``--score gate`` / :data:`mhcmatch.rank.GATE`.

.. automodule:: mhcmatch.rank
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.precursor module
-------------------------

Optional extra: ``pip install 'mhcmatch[precursor]'``.

Since 0.12.0 this is a **re-export of** ``vdjmatch.precursor`` — the estimators, their maths and the
``vdjmatch precursor`` CLI live in the repertoire library, which is where that half of the problem
belongs. The name is kept so existing imports and notebooks keep working, and
``from mhcmatch import precursor as P`` still gives ``P.event_ratio``, ``P.observed_mass``,
``P.coverage_corrected_mass``, ``P.ball_mass``, ``P.shell_profile``, ``P.motif_mass`` and
``P.cross_check``, plus what vdjmatch added: ``union_mass``, ``closed_ball_mass``,
``unseen_junctions`` and ``precursor_frequency``.

There is deliberately no ``automodule`` here: the module has no API of its own, so autodoc would
either duplicate vdjmatch's reference or, without the extra installed, fail the build. See
vdjmatch's own documentation for the signatures.

mhcmatch.data.aa_tables module
------------------------------

The vendored amino-acid property tables. Their values are not reproduced here — they are plain
``dict[str, float]`` and large; :doc:`property_basis` states what their principal components are.

.. automodule:: mhcmatch.data.aa_tables
   :members:

mhcmatch.data.contact_profile module
------------------------------------

.. automodule:: mhcmatch.data.contact_profile
