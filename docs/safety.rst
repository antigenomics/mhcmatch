Safety, prior evidence, and what goes in the cassette
=====================================================

Ranking says which mutations are worth targeting. This page is the step after: **what to withdraw
on safety grounds, what prior evidence a candidate already carries, and how many units to spend.**
They are three different questions and none of them is answered by the candidate's score.

Everything here is computed, not looked up. That matters for the one caveat this page repeats:
on a corpus assembled *from* the reference sets, the prior-evidence columns are self-fulfilling and
tell you nothing; on a fresh patient's variants they are the most informative columns in the table.

.. contents::
   :local:
   :depth: 1


The safety screen: two clauses, and it excludes rather than down-ranks
----------------------------------------------------------------------

:func:`mhcmatch.vector.screen` withdraws a unit outright. The second-best cassette is cheap and
myocarditis is not, and capacity spent on a unit that has to be withdrawn is capacity not spent on
a safe one.

The precedent is not hypothetical and it was not a binding-prediction failure:

* An affinity-enhanced TCR against the HLA-A\*01:01 MAGE-A3 epitope ``EVDPIGHLY`` killed the first
  two patients infused, by cardiogenic shock within days. Autopsy found myocardial damage with **no
  MAGE-A3 expressed in heart at all**; the off-target was ``ESDPIVAQY``, from titin (Linette *et
  al.*, *Blood* 2013;122(6):863-71, PMID 23770775; Cameron *et al.*, *Sci Transl Med*
  2013;5(197):197ra103, PMID 23926201).
* A TCR recognising MAGE-A3/A9/A12 caused necrotising leukoencephalopathy and two deaths, because
  MAGE-A12 turned out to be transcribed in human brain (Morgan *et al.*, *J Immunother*
  2013;36(2):133-51, PMID 23377668).

Both were invisible to binding prediction and visible in **expression**. So the check joins a
peptide to a *protein* to a *tissue*, and it asks two questions:

``target gene``
   Is the unit's own gene transcribed in a tissue that must not be attacked? The MAGE-A12 case. No
   register search is needed to see it — the floor is 0.25 TPM rather than the conventional 5,
   because MAGE-A12 in brain caudate is **0.33**.

``unrelated self origin``
   Does any register of the unit coincide with a self peptide from a **different** gene, and is that
   gene transcribed in an essential tissue? The titin-shaped case.

**The unit's own gene is excluded from the second clause, and it has to be.** A vaccine unit is a
long window of native context: its flanking registers *are* self peptides from its own parent
protein, and its mutated register sits one substitution from that protein's wild type. Screened
naively, every unit of every cassette fires. Those matches are also the ones tolerance already
covers — the flanks are presented in normal tissue daily.


Near-exact identity, not similarity
-----------------------------------

The obvious alternative — score each register with :mod:`mhcmatch.mimicry` and flag whatever
resembles a tolerance-side reference — was built and measured against this one, on 1,000 viral
epitopes (which cannot be self, so every firing is a false positive) and 1,000 thymic peptides from
essential-tissue genes:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20

   * - route
     - false pos.
     - true pos.
   * - mimicry, anchor-masked
     - 0.693
     - 0.944
   * - self origin (this one)
     - **0.020**
     - 0.940

**Equal sensitivity, 35× the false positives.** Anchor-channel similarity to a presented reference
is *presentation*, not recognition, so a masked match fires for every peptide sharing the allele's
motif — the influenza epitope ``GILGFVFTL`` draws 14 essential-tissue hits. A screen that flags
two-thirds of a candidate list excludes nothing in practice, because nobody withdraws two-thirds of
a cassette. ``bench/results/vector_safety_screen.md``.


Why the search radius is zero
-----------------------------

``max_subs=0`` — exact coincidence — because **the decision is per unit while the search is per
register, and that multiplies.** A 27-mer carries ~70 class-I registers and is withdrawn if any one
fires, so a per-register false-positive rate that reads as small is not the rate a cassette
experiences. Measured on six hazard-free 27-mers plus one burying the real titin epitope
(``bench/results/vector_screen_radius.md``) — **units falsely withdrawn, out of six**:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - ``max_subs``
     - 9-mers
     - 9-10-11
     - 8-9-10-11
   * - 0
     - 0
     - 0
     - 0
   * - 1
     - 1
     - 1
     - **4**

Radius 0 is clean at every length set; radius 1 is clean at none. It collapses once 8-mers enter: an
8-mer plus its 152 one-substitution neighbours is ~153 of 20\ :sup:`8` sequences against the
proteome's ~68 M windows, so a chance hit per register is expected and the ~20 8-mer registers in a
unit make it near-certain. **Every setting still catches the titin unit**, so radius 1 buys nothing
and costs most of the cassette.

.. warning::

   **What this does not catch**, stated because a safety screen that oversells itself is worse than
   none. It would not have caught the titin event *as it happened*: the construct contained MAGE-A3,
   whose profile is clean (13.4 TPM testis, 0.00 elsewhere), the cross-reactive titin peptide was
   never in it, and four TCR-facing substitutions separate the two — no distance threshold reaches
   one from the other. The affinity-enhanced TCR was the cause. What this catches is the adjacent
   and commoner failure.


.. note::

   **This screen is class I, and deliberately.** The two clauses are built on a CD8 mechanism: a
   register that *is* an essential-tissue self peptide, and a target gene transcribed where it was
   assumed silent. The class-II analogue is not the same question rewritten for longer peptides ---
   CD4 self-reactivity runs through help, hypersensitivity and allergy rather than through direct
   cytotoxicity, and a screen that reused these thresholds on class-II ligands would be asserting an
   equivalence nobody has measured. ``mhcmatch vector`` accepts ``--cls mhc2`` for the *register
   vocabulary* used in junction scanning; it does not claim the exclusion policy transfers, and the
   shipped subworkflow builds a cassette from class I only.


Prior evidence: near-exact matches to known antigens
----------------------------------------------------

A candidate that is already a confirmed neoantigen, or sits one substitution from one, carries
direct evidence no model output can match. :mod:`mhcmatch.known` reports **exact** membership of
five reference sets and :mod:`mhcmatch.rank` floats those candidates into a tier of their own, with
the model score still shown beside them — burying "this is a confirmed NCI neoantigen" inside a
weighted sum lets a mediocre score dilute the one piece of direct evidence in the row.

``mhcmatch neoag`` and ``mhcmatch rank --annotate`` widen that to **near-exact**, reporting
``neoag_distance`` (substitutions to the nearest tested neoantigen), ``neoag_nearest`` (which one)
and ``neoag_n_within``. Fuzzy beats exact: held out, matching at ≤2 substitutions roughly doubles to
triples the recall of a fresh cohort's true positives over exact lookup, which is why
``--max-subs`` defaults to 2 and 1 is the tighter reading.

.. warning::

   **These columns are meaningless on a benchmark corpus and informative on a patient.** The
   reference sets *are* the deposits our own evaluation arms are built from, so a candidate drawn
   from those arms is at distance 0 to itself and the column is all-positive by construction. It
   only starts carrying information when the candidates come from somewhere the reference has never
   seen — which is exactly the delivery case. Read a filled ``neoag_distance`` column on internal
   data as a schema check, never as a result.

The same holds for the ``self`` and ``thymus`` mimicry channels, for the opposite reason: a hit
argues *reduced* immunogenicity (reactive T cells were plausibly deleted during negative selection)
**and** cross-reactivity risk, and those two conclusions pull in different directions. They are
reported, never folded into the score.


Coverage: how many units, and on which allotypes
-------------------------------------------------

:func:`mhcmatch.vector.select` spends capacity with a **per-allotype stopping rule**, not a global
budget. The clinical numbers (20 for autogene cevumeran, 34 for mRNA-4157, 20 in four pools for
NeoVax) derive from no published objective function, but the *shape* of the competition is
established: it is for the antigen-presenting cell rather than for MHC, it is strongest within an
allotype, and it can be net-positive across allotypes.

So the expected yield is a sum of independently saturating per-allotype terms, and the rule is one
line — **keep adding to an allotype while the next candidate's probability beats that allotype's
current expected yield per slot**:

.. math::

   E = \sum_a \frac{n_0\,S_a(n_a)}{n_0 + n_a},
   \qquad
   \text{add to } a \iff p_{a,n+1} > \frac{S_a(n_a)}{n_0 + n_a}

Because each allotype saturates on its own, a crowded allotype's next unit falls below an empty
allotype's first one, so **diversification across the patient's allotypes falls out of the
arithmetic** instead of being imposed as a quota. ``n0`` is per-allotype capacity, it is the one
free parameter, and it has **no default** — nothing in the public record fits it, so the value is
the caller's to defend and :attr:`~mhcmatch.vector.Selection.n0` carries it into the result so a
cassette can always name its own assumption.

Screening runs **before** selection, for the same reason: capacity spent on a unit that will be
withdrawn is capacity not spent on a safe one.


Running it
----------

.. code-block:: bash

   mhcmatch vector \
       --candidates ranked.tsv \
       --context windows.fasta \
       --n0 6 --alleles A*02:01,B*07:02 --cls mhc1 \
       --screen \
       --out cassette.tsv --fasta cassette.faa --fasta-nt cassette.fna

``--candidates`` and ``--context`` are both required and neither is redundant: ``rank`` emits
**minimal epitopes** and a vaccine unit is the long window around the mutation. Injecting a minimal
epitope is not a smaller version of the right thing — a 9-mer loads onto any cell without
costimulation and is the tolerising configuration (PMID 17911588) — so the reader refuses a table it
cannot tell apart rather than guessing.

``cassette.tsv`` carries a ``withdrawn`` section, one row per (unit, register, source gene) with the
clause that fired, the gene, the substitutions, the tissue and its TPM. **A withdrawn candidate has
to say what withdrew it**: "the screen dropped 3 of 40" is not a safety argument, and the reason is
what a clinician overrides or accepts.

.. note::

   ``--screen`` builds one whole-proteome index **per register length** — ~12 GB peak each and a few
   minutes apiece — so screen the whole candidate list in one process, never one unit per
   invocation. **Without the flag no safety check runs at all** and the cassette carries whatever it
   was handed. On a cluster this is the process that needs 48 GB; see
   ``integrations/nextflow/mhcmatch/README.md``.

A site with its own toxicity list substitutes the policy wholesale: ``screen`` takes
``risk(unit, registers) -> [reason, ...]`` as an argument, so the shipped
:func:`~mhcmatch.vector.self_origin_risk` is one implementation and not the interface.
