# Make `modelName` report the model/facade, not the contrast-test schema

## Motivation

A user (WU347682, `prolfquapp.DIANN`, `-m rfit`) ran the **rfit** model and was
confused that the 2-group report / Excel showed:

| facade | modelName          |
|--------|--------------------|
| rfit   | WaldTest_moderated |

The chosen backend (`rfit`) is only visible in the `facade` column. The
`modelName` column reports `WaldTest_moderated` — the **contrast-test schema**,
which is identical for `lm`, `rlm`, `rfit`, `firth`, ... because they all route
through `Contrasts -> ContrastsModerated`. So `modelName` carries no information
about *which model was actually fit*, and reads as if a different model was used
than the one requested on the command line.

## Goal

1. **`modelName` should name the model / facade actually used**: `lm`,
   `lm_impute`, `rlm`, `rfit`, `firth`, `limma`, ... — i.e. mirror the `facade`
   column's information (or replace the need for a separate `facade` column).
2. **The "moderated Wald test" wording moves to the report _text_** (methods
   section), where it belongs as a description of the testing procedure shared
   by all Wald-type facades — not as a per-row column value.
3. **Document the change** in the report vignette and **add a test** that pins
   the expected `modelName` values per facade.

## Where the current value comes from (prolfqua, NOT prolfquapp)

The `modelName` column is produced inside **prolfqua**, not prolfquapp:

- `prolfqua/R/Contrasts.R:121` — default `model_name = "WaldTest"`; written to
  the `modelName` column at `Contrasts.R:225`.
- `prolfqua/R/ContrastsModerated.R:70` — wraps and appends `_moderated`
  -> `"WaldTest_moderated"`.
- `prolfqua/R/ContrastsModerated.R:138-150` — imputed/rescued rows are retagged
  by string surgery: `*_imputed` -> `*_moderated_imputed`, else
  append `_moderated`. So **`modelName` currently encodes TWO things at once**:
  the test schema AND whether the row was imputed.
- The facade key (`lm`, `rfit`, ...) is added separately as the **`facade`**
  column by `.add_facade_column()` (`prolfqua/R/ContrastsFacades.R:45`), called
  per facade (e.g. rfit at `:656`, lm at `:517`, rlm at `:585`).

Because the value originates in prolfqua, decide the fix location:

- **Option A (prolfqua, preferred):** make each facade pass a meaningful
  `model_name` into `Contrasts$new(..., model_name = <facade>)` so `modelName`
  becomes the facade name natively, and drop the `_moderated` suffix logic from
  the column (keep it as a config/methods string). Cleanest, benefits all
  downstream consumers (prophosqua, prolfquasaint, benchmark). Risk: changes a
  long-standing public column value — coordinate across the ecosystem.
- **Option B (prolfquapp-only):** after `get_contrasts()`, overwrite
  `modelName <- facade` in the output-assembly step in prolfquapp. Lower blast
  radius, but leaves prolfqua inconsistent and duplicates the `facade`/`modelName`
  columns.

Recommend **Option A** with prolfquapp follow-ups below.

## ⚠️ Do not lose the imputation flag

This is the crux. Today `modelName` doubles as the "was this row imputed?"
signal (`WaldTest_moderated_imputed` / `Imputed_Mean_moderated`). Several
consumers rely on it:

- **Volcano coloring**: `prolfqua/R/ContrastsPlotter.R:171` defaults the point
  `colour` to the `modelName` column — imputed proteins are drawn as gray dots.
- **Color palette**: `prolfquapp/R/report_helpers.R:150` `.model_name_palette()`
  hard-codes keys `Linear_Model_moderated`, `Imputed_Mean_moderated`,
  `WaldTest_moderated`. Changing `modelName` values **breaks these defaults** —
  must be updated.

If `modelName` becomes purely the facade name, **add a separate column** (e.g.
`imputed` boolean, or keep an `estimate_type` column) so the observed-vs-imputed
distinction survives, and repoint the volcano `colour` and the palette at it.

## Concrete change list

### prolfqua (if Option A)
- [ ] `R/ContrastsFacades.R` — each facade passes its key as `model_name` to
      `Contrasts$new()` (or via `ContrastsModerated`). Confirm: `lm`, `lm_impute`,
      `lm_missing`, `rlm`, `rfit`, `firth`, `limma`(+voom/impute), `deqms`,
      `limpa`, and nested facades (`lmer_nested`, `ropeca_nested`, ...).
- [ ] `R/ContrastsModerated.R:70,138-150` — stop folding `_moderated` /
      `_moderated_imputed` into `modelName`; surface imputation as its own flag.
- [ ] `R/Contrasts.R:121` — reconsider the `"WaldTest"` default.
- [ ] Update `test-ContrastsFacades.R`, `test-Contrasts.R`,
      `test-ContrastsModerated*` to assert the new `modelName` per facade.

### prolfquapp
- [ ] `R/report_helpers.R:150` `.model_name_palette()` — rekey the default
      palette to facade names (or to the new imputation flag), so report colors
      stay stable.
- [ ] `R/se_report_lfqdata.R:302-303` and
      `inst/templates/quarto/Grp2Analysis_V2_SE.qmd:223` — they map/default
      `modelName`; check they still make sense.
- [ ] Move the moderated-Wald-test explanation into the **methods text** of the
      report (already partly present at
      `vignettes/Grp2Analysis_V2_R6.Rmd:101`: "protein variances are
      moderated [@Smyth2004]...").

### Vignette / docs
- [ ] `vignettes/Grp2Analysis_V2_R6.Rmd:359,373` currently describe `modelName`
      as "Imputed_mean or Linear_Model_Moderated" — **already stale** (the actual
      value is `WaldTest_moderated`). Rewrite to: `modelName` = the model used
      (`lm`, `rfit`, ...), plus a sentence in the methods section stating all
      these backends are tested with an empirical-Bayes **moderated Wald test**
      (Smyth 2004), and that imputed proteins are flagged separately and drawn
      in gray.
- [ ] Add/extend a test (prolfquapp `tests/testthat/`) that runs a small DEA per
      facade and asserts the `modelName` column equals the requested facade key.

## Acceptance criteria

- For `-m rfit`, the report/Excel `modelName` column reads `rfit` (not
  `WaldTest_moderated`); same pattern for `lm`, `lm_impute`, `rlm`, ...
- Imputed proteins are still distinguishable (separate flag) and still render as
  gray dots in the volcano; report colors unchanged or intentionally updated.
- The "moderated Wald test (Smyth 2004)" description appears in the report
  methods text.
- Vignette `Grp2Analysis_V2_R6.Rmd` text matches the new column semantics.
- Tests pin `modelName` per facade and pass.

## Open questions

- Option A vs B (change prolfqua vs override in prolfquapp)?
- Keep the `facade` column at all if `modelName` now carries the same info, or
  drop one to avoid redundancy?
- Name of the new imputation flag column (`imputed` boolean vs a categorical
  `estimate_type` of observed/imputed)?
