# AGENTS.md

Guidance for AI coding agents working in this repository.

- Read `docs/DESIGN.md` before changing behaviour. Paradise is not a physics solver: people
  choose the plate motions, and the engine does the bookkeeping. Do not add forces that choose
  motions or events on the person's behalf. The advisor in the plan may propose them, but never
  apply them.
- `paradise/core` is copied from Poly-Time. Record every change to it in
  `paradise/core/PROVENANCE.md`, and keep its inherited tests (`tests/test_core_*.py`) passing.
- Sub-steps are a fixed 2.5 Myr by the developer's decision. Do not shorten them for accuracy.
- Run the tests with `python -m pytest -q tests` (about two minutes).
