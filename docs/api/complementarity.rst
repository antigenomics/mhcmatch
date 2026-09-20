Complementarity --- the I and C of EPIC
---------------------------------------

The recognition axis: chemistry, identity and the corpus channels.

mhcmatch.complement module
~~~~~~~~~~~~~~~~~~~~~~~~~~

The recognition axis as one score: whole-peptide physicochemistry and length, the same components
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The head dispatcher over the recognition axis: ``complement`` (the six-block model above, and the
default), plus ``posbayes``, ``physchem_glm`` and ``esm64_glm``, each fitted alone so their BIC is
comparable. See :doc:`complementarity`.

.. automodule:: mhcmatch.recognition
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.immuno module
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mhcmatch.immuno
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.posbayes module
~~~~~~~~~~~~~~~~~~~~~~~~

Position-role naive Bayes over amino-acid identity: anchor and TCR-facing residues get separate
conditional distributions, because for several amino acids their contributions carry opposite signs
and pooling averages that away. Emits a **log-likelihood ratio**, so the caller supplies the prior.

Grouped 5-fold CV **0.712** human / **0.758** mouse;
size-matched transfer **0.731** human→mouse, **0.692** mouse→human. Cysteine is masked — see the
module warning.

.. automodule:: mhcmatch.posbayes
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.luksza module
~~~~~~~~~~~~~~~~~~~~~~

The Łuksza recognition term :math:`R = Z/(1+Z)` -- a soft partition function over near-matches,
replacing a hard distance cut, so **how many** near-matches a candidate has and **how near** they
are both enter. ``viral_R`` is not a term of the shipped aggregate — the sum
lived only in the benchmark repo, which made :func:`mhcmatch.rank.aggregate_score` callable with a
feature no installed user could supply. ``EPIC`` does not score with it; the quantity is still the
published recognition term and is
still computed here. ``k`` and ``a0`` are read from the shipped artifact when one carries them, so
a refit needs no code change.

.. automodule:: mhcmatch.luksza
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.mimicry module
~~~~~~~~~~~~~~~~~~~~~~~

Mimicry as a signed, per-component immune-response risk. Three references — ``viral`` (priming),
``self`` (tolerance, and the autoimmunity read-out) and ``thymus`` (negative selection) — each split
into an **anchor** and a **TCR-facing** channel, because a whole-peptide distance averages two
different measurements. Scores are **log-odds**; :func:`mhcmatch.mimicry.probability` is a separate
step that requires a *named* corpus, because the screens behind any calibration run from 0.048 % to
46.8 % positive and an unqualified probability is mostly a statement about which intercept was used.

The tested-neoantigen database is exposed as :func:`mhcmatch.mimicry.annotate` — **prior evidence,
never a fitted term**. Every labelled screen we hold sits inside that database, so a coefficient on
it would be memorisation; held out of the fit, fuzzy matching at two substitutions still recovers
0.08–0.34 of a fresh screen's positives against 0.00–0.26 for exact lookup, which is what makes it
worth reporting. Command line: ``mhcmatch neoag``.

.. automodule:: mhcmatch.mimicry
   :members:
   :undoc-members:
   :show-inheritance:

mhcmatch.mimics module
~~~~~~~~~~~~~~~~~~~~~~

Cross-reactivity by **category**, never summed: a hit in the thymic immunopeptidome argues tolerance
and autoimmune risk, a hit in the host proteome argues peripheral cross-reactivity without implying
presentation, and a viral or bacterial hit argues a pre-existing repertoire may cross-react.
:data:`~mhcmatch.mimics.KINDS` states which is which; :func:`~mhcmatch.mimics.neighbours` runs the
whole query set through one threaded C++ index per (category, length).

.. automodule:: mhcmatch.mimics
   :members:
   :undoc-members:
   :show-inheritance:
