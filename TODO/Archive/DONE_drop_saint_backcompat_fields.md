# TODO: drop transitional `dea$saint_input` / `dea$saint_result` fields

Created: 2026-06-16

Carved out of the now-completed SAINT facade integration
(`TODO/Archive/DONE_saint_facde_integration.md`, Phase D3) as the one
deliverable intentionally deferred "until after one release".

## Background

SAINT is now reached generically through `prolfqua::lookup_facade("saint")`,
and report code reads SAINT artifacts via `contrast_obj$extra_artifacts()`
(`R/R6_DEAReportGenerator.R`, `prep_result_list()` merges that list).

For one-release back-compatibility, `DEAnalyse` still mirrors the facade's
extras onto two transitional fields so any external consumer reading
`dea$saint_input` / `dea$saint_result` keeps working:

- `R/R6_DEAnalyse.R:84-87` — the `saint_input` / `saint_result` `@field`s.
- `R/R6_DEAnalyse.R:177-190` — the block that copies `facade$extra_artifacts()`
  (`saint_inter`/`saint_prey`/`saint_bait`/`saint_list`) onto `self$saint_input`
  / `self$saint_result` after `build_facade()`.

## What to do (when the deferral window has passed)

- Remove the `saint_input` / `saint_result` fields and the mirroring block in
  `R/R6_DEAnalyse.R`.
- Confirm nothing internal still reads `dea$saint_input` / `dea$saint_result`
  (grep `R/`, `inst/`, `vignettes/`, `tests/`); everything should go through
  `contrast_obj$extra_artifacts()`.
- `make document` to regenerate the `.Rd`, then `make check-fast` and the
  saint tests (`tests/testthat/test-saint-model.R`).

## Notes / precondition

- The original plan said to **deprecate the fields in roxygen first**, then
  remove them after a release. The current `@field` docs are plain (no
  lifecycle/deprecation note), so either add a deprecation note now and remove
  next cycle, or confirm there are no external consumers and remove directly.
- This is cleanup only — no behavior change for the supported
  `extra_artifacts()` path.
