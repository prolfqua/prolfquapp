# Merge `run_dea_peptide_to_protein` into `run_dea`

## Interface (live in prolfqua)

`prolfqua::lookup_facade(name)$needs` returns exactly two values:

- `"same"` — facade emits contrasts at the same hierarchy level as its input
  (protein → protein FC, or peptide/precursor → peptide/precursor FC).
- `"nested"` — facade takes child-level input and emits parent-level
  contrasts (peptide/precursor → protein FC).

### Reader / facade compatibility

| Reader (software key)                                                                                   | Valid facade `needs` |
|---------------------------------------------------------------------------------------------------------|----------------------|
| `DIANN`, `MAXQUANT`, `FP_TMT`, … (protein-aggregated)                                                   | `"same"` only        |
| `DIANN_PEPTIDE`, `FP_TMT_PEPTIDE`, `MAXQUANT_PEPTIDE`, `MSSTATS_PEPTIDE`, `MSSTATS_FP_DIA_PEPTIDE`, `BGS_PEPTIDE` | `"same"` or `"nested"` |

Enforced at the top of `run_dea`.

## Scope (decided)

- Single registry; dispatch on `lookup_facade(model)$needs`.
- Delete `run_dea_peptide_to_protein()`, `CMD_DEA_PEPTIDE_TO_PROTEIN.R`, and the
  shell wrapper outright. No deprecation period.
- Full sweep: function + .Rd + CMD R script + shell wrapper + tests + NAMESPACE.

## Changes

### 1. `prolfquapp/R/cmd_helpers.R`

Merge the two functions into a single `run_dea(indir, dataset, software, config)`:

- Resolve the facade and its registry entry once at the top:
  ```r
  default_model <- .resolve_facade_model(
    config$processing_options$model,
    config$processing_options$model_missing
  )
  model_entry <- prolfqua::lookup_facade(default_model)
  if (is.null(model_entry)) stop("Unknown facade: ", default_model)
  is_nested <- identical(model_entry$needs, "nested")
  saint_annot <- isTRUE(model_entry$needs_saint_annotation)
  ```
- Enforce reader / facade pairing:
  ```r
  is_peptide_reader <- grepl("_PEPTIDE$", software)
  if (is_nested && !is_peptide_reader) {
    stop("Facade '", default_model, "' (needs=\"nested\") requires a ",
         "peptide-level reader (e.g. DIANN_PEPTIDE). Got '", software, "'.",
         call. = FALSE)
  }
  ```
  Protein reader + `"same"` facade is the legacy default. Peptide reader +
  `"same"` facade is also legal — yields peptide-level fold changes.
- Read annotation with `SAINT = saint_annot` (current `run_dea` behaviour;
  resolves to `FALSE` for non-SAINT facades).
- Validate `software` against `get_procfuncs()` and run
  `preprocess_software(...)` (unchanged).
- `data_prep <- ProteinDataPrep$new(xd$lfqdata, xd$protein_annotation, config)`,
  then `cont_decoy_summary()` and `remove_cont_decoy()` (shared).
- Branch:
  - **`is_nested`** — port the peptide-keep path from current
    `run_dea_peptide_to_protein` (the `"nested"` sub-branch of lines 344–396 in
    `cmd_helpers.R`). The `"aggregated_limpa"` and `"either"` sub-branches and
    the trailing `stop()` fallback all become dead code under the two-value
    interface — delete them. What remains:
    ```r
    lfq_peptide <- data_prep$lfq_data_peptide$get_copy()
    lfq_peptide$set_config_value("hierarchy_depth", 1)
    lfq_raw <- lfq_peptide
    lfq_model <- prolfquapp::transform_lfqdata(
      lfq_peptide,
      method = config$processing_options$transform
    )
    report_prep <- prolfquapp::ProteinDataPrep$new(
      lfq_peptide$get_copy(), data_prep$rowAnnot, config
    )
    report_prep$aggregate()
    report_prep$transform_data()
    lfq_raw$rename_response("abundance")
    lfq_model$rename_response("normalized_abundance")
    deanalyse <- DEAnalysePeptideToProtein$new(
      lfq_data = lfq_model,
      rowAnnot = data_prep$rowAnnot,
      prolfq_app_config = config,
      contrasts = annotation$contrasts,
      default_model = default_model,
      lfq_data_raw = lfq_raw,
      summary = data_prep$summary
    )
    deanalyse$build_default()
    deanalyse$get_annotated_contrasts()
    deanalyse$lfq_data <- report_prep$lfq_data_transformed
    deanalyse$lfq_data_raw <- report_prep$lfq_data
    ```
  - **`!is_nested`** (i.e. `"same"`) — existing protein path:
    `data_prep$aggregate(); data_prep$transform_data();
    data_prep$build_deanalyse(annotation$contrasts);
    deanalyse$build_default(); deanalyse$get_annotated_contrasts()`.
- Return `list(deanalyse, xd, annotation, files)` in both branches.

Delete the body of `run_dea_peptide_to_protein()` (lines 284–404).

### 2. `prolfquapp/NAMESPACE`

`make document` removes `export(run_dea_peptide_to_protein)` automatically
once the source function is gone.

### 3. `prolfquapp/man/run_dea_peptide_to_protein.Rd`

Deleted by `devtools::document()`. Verify after `make document`.

### 4. `prolfquapp/inst/application/CMD_DEA_PEPTIDE_TO_PROTEIN.R`

Delete. `CMD_DEA.R` handles both `"same"` and `"nested"` facades since
dispatch lives in `run_dea`.

Sanity-check that `CMD_DEA.R`:

- Accepts the `_PEPTIDE` software keys (already listed by `get_procfuncs()`).
- Accepts nested facade names (`limpa_nested`, `lmer_nested`, `ropeca_nested`,
  `firth_nested`) via YAML `processing_options$model` or the `-m/--model` flag
  through `sync_opt_config()`.

### 5. `prolfquapp/inst/application/bin/prolfqua_dea_peptide_to_protein.sh`

Delete the shell wrapper.

### 6. `copy_shell_script()` registry

```
grep -rn "prolfqua_dea_peptide_to_protein\|dea_peptide_to_protein" prolfquapp/R/
```

Remove the entry so `copy_shell_script(workdir = ".")` doesn't try to copy a
deleted file.

### 7. Tests

```
grep -rn "run_dea_peptide_to_protein\|CMD_DEA_PEPTIDE_TO_PROTEIN\|prolfqua_dea_peptide_to_protein" \
  prolfquapp/tests/ prolfquapp/inst/
```

- Replace `prolfquapp:::run_dea_peptide_to_protein(...)` with
  `prolfquapp::run_dea(...)`.
- If any test invokes `CMD_DEA_PEPTIDE_TO_PROTEIN.R` as a subprocess, switch to
  `CMD_DEA.R`.
- `DEAnalysePeptideToProtein` (R6 class) is still used internally by the
  nested branch — keep the class and its export.

### 8. Documentation

```
grep -rn "peptide_to_protein\|PEPTIDE_TO_PROTEIN\|aggregated_limpa\|\"aggregated\"" \
  prolfquapp/vignettes/ prolfquapp/AGENTS.md prolfquapp/README* 2>/dev/null
```

- `prolfquapp/AGENTS.md` "Command-Line Scripts" section: remove any bullet for
  `CMD_DEA_PEPTIDE_TO_PROTEIN.R`.
- Vignettes mentioning the peptide-to-protein script: redirect to `CMD_DEA.R`
  with a nested facade in the YAML.
- Replace any lingering `"aggregated"` / `"aggregated_limpa"` references with
  the two-value contract (`"same"` / `"nested"`).

## Verification

1. `cd prolfquapp && make document` — confirms
   `prolfquapp/man/run_dea_peptide_to_protein.Rd` is gone and `NAMESPACE` no
   longer lists `run_dea_peptide_to_protein`.
2. `make install` — package builds.
3. `make test` — updated tests pass.
4. `make check-fast` — no orphaned `.Rd` / NAMESPACE entries.
5. Manual smoke tests:
   - Protein reader + same facade: `CMD_DEA.R --software DIANN`, YAML
     `model: limma`. Should match the pre-refactor result.
   - Peptide reader + nested facade: `CMD_DEA.R --software DIANN_PEPTIDE`,
     YAML `model: limpa_nested`. Should match the pre-refactor result from
     the deleted peptide-to-protein path.
   - Peptide reader + same facade: `CMD_DEA.R --software DIANN_PEPTIDE`,
     YAML `model: limma`. Peptide-level fold changes (new path).
   - Protein reader + nested facade: `CMD_DEA.R --software DIANN`, YAML
     `model: limpa_nested`. Must error with the pairing message.

## Out of scope

- Renaming `_PEPTIDE` software keys.
- Auto-promoting a protein reader to its peptide variant when a nested facade
  is selected.
- Further prolfqua changes — the two-value `needs` interface is already in
  place.
