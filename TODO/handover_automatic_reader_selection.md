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

## Open question — is the slurmworker a better place?

Quite possibly **yes**, or at least *also*. Two layers could own this:

1. **prolfquapp `run_dea()` (done here).** Pro: the invariant "nested facade ⇒
   peptide reader" lives next to the code that knows about facades and readers;
   every caller (CLI, tests, ad-hoc) benefits; it is unit-testable. Con: it
   silently changes the effective reader at run time; the workunit definition on
   B-Fabric still records the "wrong" `software`, so provenance/params don't
   reflect what actually ran.

2. **slurmworker / app config** (`/home/bfabric/slurmworker/config/A414_DEA/app.yml`
   and the dispatch layer that builds `opt$software` for `CMD_DEA_V2.R`). Pro:
   the correct `software` could be written into the workunit's resolved
   parameters *before* the job runs, so the recorded params match execution
   (better provenance); it could also surface the choice back to B-Fabric. Con:
   the slurmworker would need to know the facade→`needs` mapping (i.e. call
   `prolfqua::lookup_facade`) and the `_PEPTIDE` naming convention, duplicating
   knowledge that currently lives in prolfqua/prolfquapp.

**Recommendation.** Keep the `run_dea()` safety net (it guarantees correctness
regardless of how the job was launched), and *additionally* consider resolving
the reader at the slurmworker/dispatch layer so the **recorded** workunit
parameters match what runs. The cleanest split:

- prolfquapp keeps the authoritative mapping logic (exporting a tiny helper the
  worker can call, e.g. a public wrapper around `.resolve_nested_reader`, so the
  convention is defined in exactly one place);
- the slurmworker calls that helper during dispatch to normalise `software`
  *and* writes the normalised value into the workunit params/log.

That avoids duplicating the `_PEPTIDE` convention in two repos while fixing the
provenance gap that the in-`run_dea` approach leaves open.

### Follow-up tasks

- [ ] Decide whether to expose `.resolve_nested_reader` as an exported function
      for the slurmworker to reuse.
- [ ] If yes: implement reader normalisation in the slurmworker dispatch and
      persist the resolved `software` into the workunit definition/params.
- [ ] Consider whether `prolfqua_dataset.sh` / YAML generation should also warn
      when a nested model is configured against a protein-level reader, so users
      see the implication at configuration time rather than only at run time.
