Optional extras
---------------

Re-exports of a sibling library, gated behind an install extra.

mhcmatch.precursor module
~~~~~~~~~~~~~~~~~~~~~~~~~

Optional extra: ``pip install 'mhcmatch[precursor]'``.

This is a **re-export of** ``vdjmatch.precursor`` — the estimators, their maths and the
``vdjmatch precursor`` CLI live in the repertoire library, which is where that half of the problem
belongs. The name is kept so existing imports and notebooks keep working, and
``from mhcmatch import precursor as P`` still gives ``P.event_ratio``, ``P.observed_mass``,
``P.coverage_corrected_mass``, ``P.ball_mass``, ``P.shell_profile``, ``P.motif_mass`` and
``P.cross_check``, plus what vdjmatch added: ``union_mass``, ``closed_ball_mass``,
``unseen_junctions`` and ``precursor_frequency``.

There is deliberately no ``automodule`` here: the module has no API of its own, so autodoc would
either duplicate vdjmatch's reference or, without the extra installed, fail the build. See
vdjmatch's own documentation for the signatures.
