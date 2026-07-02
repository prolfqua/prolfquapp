# Review: expose `nr_peptides` end-to-end

Date: 2026-07-02
Reviewed plan: `TODO_expose_nr_peptides.md` after the original review was folded in.

This is a static review against the current working tree, plus one small empirical check of the existing
`prolfqua::LFQData$filter_proteins_by_peptide_count()` behavior. I did not run package tests.

## Summary

The revised plan fixes the main issues from the previous review. It no longer proposes passing `nr_peptides` through
every `do.call()` adapter, it records A414 YAML emission as done, and it correctly keeps `ProteinAnnotation` out of the
filtering decision.

The remaining risks are now about the exact insertion point and helper contract for the generic LFQData-side filter.
The plan should be tightened before implementation so the filtered LFQData, returned `xd`, `ProteinAnnotation`, IBAQ,
and reports cannot drift apart.

## Findings

### 1. High: filtered LFQData and `ProteinAnnotation` can become inconsistent unless the re-attachment step is specified

The plan says to call the generic helper from `run_dea()` and then "re-attach `ProteinAnnotation` (right-join)
afterward" (`TODO_expose_nr_peptides.md`:98-100). That is the right concept, but it is not yet precise enough for the
current code path.

Today every reader returns both `xd$lfqdata` and an already-constructed `xd$protein_annotation` before `run_dea()` gets
control. `run_dea()` then constructs `ProteinDataPrep$new(xd$lfqdata, xd$protein_annotation, config)`
(`R/cmd_helpers.R:330-343`). Output writing later still reads from the returned `xd`: it logs
`xd$protein_annotation$get_summary()` and computes IBAQ from `xd$lfqdata` plus `xd$protein_annotation`
(`R/cmd_helpers.R:418-472`; duplicated in `inst/application/CMD_DEA_V2.R:184-232`).

If the new helper only filters `data_prep$lfq_data_peptide` or only the `DEAnalyse` inputs, the model can be filtered
while `result$xd`, the logged protein summary, and IBAQ still reflect the unfiltered reader output. If it filters
`xd$lfqdata` but reuses the old `xd$protein_annotation`, annotation remains built against proteins that no longer
exist in the quant layer. The right-join invariant prevents row loss, but it does not make stale summaries or stale
child-count metadata correct.

Recommendation: make the implementation item explicit: after filtering, update the canonical `xd` object used by
`run_dea()` and returned to output writing. Either rebuild `ProteinAnnotation` against the filtered `LFQData` using the
old annotation table as metadata, or add a narrow resync helper that preserves `description`, `cleaned_ids`, `full_id`,
`exp_nr_children`, and patterns while rebuilding from the filtered quant rows. Add a regression asserting that
`result$xd$lfqdata`, `result$xd$protein_annotation$row_annot`, `deanalyse$lfq_data_raw`, and IBAQ all agree on the kept
proteins.

### 2. High: the existing prolfqua peptide-count filter cannot be reused unchanged

There is already an exported `prolfqua::filter_proteins_by_peptide_count()` and an `LFQData$filter_proteins_by_peptide_count()`
method (`prolfqua/R/LFQData.R:181-184`, `prolfqua/R/LFQData.R:549-560`). The TODO's open helper-location question
(`TODO_expose_nr_peptides.md`:137-140) should discuss this existing API explicitly.

The current prolfqua helper is not equivalent to the revised plan. It reads the threshold from
`AnalysisConfiguration$min_peptides_protein`, not `ProcessingOptions$nr_peptides`, and its count is based on
`config$hierarchy_keys_depth()` plus the next hierarchy level (`prolfqua/R/LFQData.R:520-528`). That means it filters a
protein-depth LFQData, but becomes a no-op for `_PEPTIDE` readers where `hierarchy_depth = 2` and there is no "next"
hierarchy level. I confirmed this empirically: with the simulated peptide LFQData, depth 1 reduced rows at threshold 2,
while depth 2 left the row count unchanged.

Recommendation: do not call the existing method blindly from `run_dea()`. Either implement a prolfquapp helper that
always counts `parent = lfq$hierarchy_keys()[1]` against an explicit stripped-peptide key, independent of modelling
`hierarchy_depth`, or refactor the prolfqua helper so it accepts explicit parent/child keys and threshold. Tests must
cover both protein-level and `_PEPTIDE`/nested readers.

### 3. Medium: unsupported downstream readers need a hard gate, not accidental counting

The plan defers PTM/site readers until they expose a stripped-peptide hierarchy element
(`TODO_expose_nr_peptides.md`:90-92, 130-132`). But `run_dea()` discovers readers dynamically through `get_procfuncs()`,
which combines preprocess function maps from installed `prolfqua*` packages (`R/find_function_in_packages.r:42-52`).
So a generic `run_dea()` filter will also see downstream readers when they are installed.

If `nr_peptides > 1` is applied to a site-level or other downstream LFQData before the stripped peptide key exists, a
generic "first hierarchy child" implementation could count sites per protein instead of stripped peptides per protein,
or silently no-op if the expected key is missing. Both would violate the stated contract.

Recommendation: add an explicit support check to the helper. For `nr_peptides <= 1`, return unchanged. For
`nr_peptides > 1`, require either a known supported software family or an explicit stripped-peptide key present in the
LFQData. Unsupported downstream readers should get a clear error or an intentionally documented skip, not accidental
site-count filtering.

### 4. Medium: the MZMine exemption should be explicit, not inferred from missing hierarchy depth

The TODO lists MZMine as exempt and leaves the mechanism open: explicit software skip vs. "naturally having no peptide
hierarchy key" (`TODO_expose_nr_peptides.md`:90, 141-143). Inferring the exemption from a missing child hierarchy would
also exempt any other single-level or already-aggregated LFQData, including future unsupported readers. That turns a
contract violation into a silent no-op.

Recommendation: use an explicit MZMine/MZMineannot skip in `run_dea()` or in the helper's support policy, and test that
MZMine is unchanged at `nr_peptides = 2`. For non-MZMine LFQData without a peptide hierarchy key, prefer a clear error
when `nr_peptides > 1`.

### 5. Low: validation of `nr_peptides` is not yet part of the plan

The plan defines default `N = 1`, adds native YAML, CLI, and generated-YAML surfaces, and A414 now writes
`processing_options$nr_peptides` (`TODO_expose_nr_peptides.md`:18-20, 102-104, 127-129). It does not yet say where the
value is validated. A hand-written native YAML or CLI override can supply `0`, `-1`, `1.5`, `NA`, or a non-numeric
string unless prolfquapp normalizes and checks it centrally.

Recommendation: add a small normalization/validation step for `processing_options$nr_peptides`: require a scalar whole
number `>= 1`, coerce acceptable YAML numerics cleanly, and fail with a clear error otherwise. Cover native YAML and CLI
override tests. This is especially important because the filtering helper should be allowed to assume a valid threshold.

## Test additions I would require

- Helper unit test: counts distinct stripped peptides per top-level protein independently of `hierarchy_depth`.
- Depth regression: the filter works for both normal protein-level readers and `_PEPTIDE`/nested readers.
- State consistency test: after `run_dea(nr_peptides = 2)`, `result$xd$lfqdata`, `result$xd$protein_annotation`, and
  `result$deanalyse$lfq_data_raw` agree on the kept proteins.
- Output-path test or focused unit: IBAQ/log-summary input uses the filtered `xd`, not stale reader output.
- MZMine test: threshold is explicitly ignored and features are unchanged.
- Unsupported-reader test: PTM/site-like LFQData without a stripped-peptide hierarchy errors or is explicitly skipped.
- Validation tests: `nr_peptides = 1`, `2`, and YAML numeric `2.0` pass; `0`, negative, non-integer, missing/NA, and
  non-numeric values fail clearly unless absent defaults to 1.

## Bottom line

The revised design is much better: one LFQData-side filter is the right direction for the DEA path, and it respects the
new `ProteinAnnotation` invariant. The plan still needs to define the exact filtered-state handoff and avoid reusing
the existing prolfqua filter unchanged, because that helper is depth-dependent and misses the `_PEPTIDE` case the TODO
explicitly wants to support.
