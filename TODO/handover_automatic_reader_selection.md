# Handover: automatic peptide-level reader selection for nested facades

## Context / motivation

Nested modelling facades (`firth_nested`, `lmer`, `ropeca`, …) fit models on
**peptide-level** data. In `prolfquapp` the quantification reader and the model
are configured independently:

- the reader (software key) decides `hierarchy_depth` — protein-level (`1`) vs
  peptide-level (`2`);
- the facade (model) decides whether nested/peptide data is required
  (`prolfqua::lookup_facade(model)$needs == "nested"`).

Until now, pairing a nested facade with a protein-level reader was a hard error.
This bit B-Fabric workunit **347649** (`/scratch/A414_DEA/WU347649`): the app was
configured with `model = firth_nested` (via `model_extra`) but
`software = prolfquapp.DIANN`, and the run aborted with:

> Facade 'firth_nested' (needs="nested") requires a peptide-level reader
> (e.g. DIANN_PEPTIDE). Got 'prolfquapp.DIANN'.

## What was done (this change, in `prolfquapp`)

Implemented automatic reader selection inside `run_dea()`
([R/cmd_helpers.R](../R/cmd_helpers.R)).

Key facts that make this safe and simple:

- In `prolfqua_preprocess_functions` ([R/preprocess_software.R](../R/preprocess_software.R)),
  every reader has a `<READER>_PEPTIDE` twin that is **identical except**
  `hierarchy_depth = 2` (e.g. `DIANN` → `DIANN_PEPTIDE`,
  `FP_TMT` → `FP_TMT_PEPTIDE`, `MAXQUANT` → `MAXQUANT_PEPTIDE`,
  `MSSTATS[_FP_DIA]` → `…_PEPTIDE`, `BGS` → `BGS_PEPTIDE`).
- `get_procfuncs()` prefixes keys with the package name, so the runtime keys are
  `prolfquapp.DIANN`, `prolfquapp.DIANN_PEPTIDE`, … The `_PEPTIDE$` suffix
  convention survives the prefixing, so `paste0(software, "_PEPTIDE")` is the
  correct mapping.

New internal helper `.resolve_nested_reader(software, is_nested, available, facade)`:

- non-nested facade, or already a `_PEPTIDE` reader → returns `software` unchanged;
- nested facade + protein reader → if `<software>_PEPTIDE` is registered, switches
  to it and logs the switch (`logger::log_info`); otherwise errors with
  `"... no peptide-level counterpart '<x>' is registered ..."`.

Readers **without** a peptide twin (`SIM`, `DUMMY`, `MZMINE`, `MZMINEannot`)
therefore still produce a clear, actionable error when combined with a nested
facade — they have no peptide-level mode.

### Tests

`tests/testthat/test-run_dea.R` covers the helper directly (switch, the two
unchanged cases, and the no-counterpart error) plus a `run_dea()`-level
no-counterpart error using the SIM fixture.

### NOT verified in the dev box

The dev machine had R 4.4.2 while the installed packages were built under R 4.5.0
(`data.table` failed to load: `undefined symbol: MAYBE_SHARED`), so prolfqua/
prolfquapp could not be loaded and the full suite could not run. The helper logic
was verified in isolation (base R + logger). **Before merging:**

- [ ] `make document` (generates the `.Rd` for the new internal helper; NAMESPACE
      is unaffected because the helper is `@keywords internal`, not exported)
- [ ] `make test` in an R 4.5 environment
- [ ] re-run WU347649 end-to-end (`make run-all` in `/scratch/A414_DEA/WU347649`)
      after reinstalling prolfquapp, to confirm the firth_nested run completes

## Follow-up on 2026-06-22

Completed the provenance follow-up in prolfquapp:

- `run_dea()` now returns both `software` (the resolved reader that actually ran) and `requested_software` (the caller's original key).
- `run_dea()` updates `config$software` after reader resolution, so downstream output writers persist the resolved reader in the app config.
- `write_dea_run_outputs()` and `CMD_DEA_V2.R` log and persist the resolved reader instead of continuing to use `opt$software`.
- The helper test no longer expects `logger::log_info()` to be a base R message condition.

Validation:

- `make document` completed and generated/updated the roxygen docs.
- `Rscript -e 'devtools::test(filter = "run_dea")'` passed: 18 passed, 0 failed, 0 skipped, 4 warnings.
- `Rscript -e 'devtools::test()'` passed before the final provenance assertion was added: 254 passed, 0 failed, 1 skipped, 28 warnings. The skipped test was the installed-package subprocess test in `test-CMD_DEA_CD.R`.
- Real-data integration coverage was added in `integration_test/tests/testthat/test-dea-diann-auto-peptide.R`, using the
  subsetted real DIA-NN WU345302 fixture. It requests `--software prolfquapp.DIANN` with `model: limpa_nested` and asserts
  that the CLI switches to `prolfquapp.DIANN_PEPTIDE`, writes a `SummarizedExperiment.rds`, and persists
  `software: prolfquapp.DIANN_PEPTIDE` in `minimal.yaml`.
- `make install` in `integration_test` completed, then `make test-dea-diann-auto-peptide` passed: 8 passed, 0 failed.

Still open:

- Re-run WU347649 end-to-end (`make run-all` in `/scratch/A414_DEA/WU347649`) after reinstalling prolfquapp.
- Optional: add `--print-resolved-software` only if the slurmworker needs the resolved reader before launching the job.
- Optional: decide whether dataset/YAML generation should warn earlier when a nested model is configured with a protein-level reader.

## Open question — is the slurmworker a better place?

Initially this looked attractive (resolve `software` at dispatch so the
**recorded** workunit params match what runs). On closer inspection it is *not*
the better place, for two reasons:

1. **The slurmworker is Python.** Owning this decision there means either a
   `python → R` subprocess round-trip just to obtain one string, or
   reimplementing the logic in Python.

2. **The CLI does not expose enough for Python to decide.**
   `prolfqua_dea.sh --help` lists the **reader names**
   (`names(get_procfuncs())`) and the **model names** (`FACADE_REGISTRY`), so
   Python *could* see that `prolfquapp.DIANN` and `prolfquapp.DIANN_PEPTIDE`
   both exist and infer the `_PEPTIDE` convention. But `--help` does **not**
   expose `needs == "nested"` — i.e. *which models require peptide data*.
   Without that, the worker cannot know *when* to remap; it would have to
   hardcode the list of nested facades, duplicating prolfqua knowledge that
   goes stale as soon as a new nested facade is registered.

So the authoritative "nested facade ⇒ peptide reader" logic belongs in
prolfquapp (where the facade→`needs` mapping and the reader registry both live),
which is what this change does. The worker layer's only real added value was
**provenance** (recording what actually ran), and that is better solved without
any cross-language coupling.

**Recommendation.**

- Keep the `run_dea()` resolution (done) as the single source of truth.
- Solve provenance *inside prolfquapp*: surface the resolved `software` in the
  run outputs so the results record what actually ran. `run_dea()` already logs
  the switch via `logger::log_info`; additionally return the resolved `software`
  in the `run_dea()` result list and persist it into the output params/YAML
  (see `write_dea_run_outputs()`). No Python↔R coupling needed.
- Only if the worker genuinely must know the reader *before* the job runs:
  expose it through the CLI (R owns the logic, Python consumes a string) — e.g.
  a `--print-resolved-software` flag on `prolfqua_dea.sh` that prints the
  resolved key and exits. Do **not** reimplement the nested-model list in
  Python.

### Follow-up tasks

- [x] Return the resolved `software` from `run_dea()` and persist it into the
      output params/YAML (`write_dea_run_outputs()`) for provenance.
- [ ] (Optional) Add a `--print-resolved-software` flag to `CMD_DEA_V2.R` /
      `prolfqua_dea.sh` if the slurmworker ever needs the resolved reader at
      dispatch time.
- [ ] Consider whether `prolfqua_dataset.sh` / YAML generation should also warn
      when a nested model is configured against a protein-level reader, so users
      see the implication at configuration time rather than only at run time.
