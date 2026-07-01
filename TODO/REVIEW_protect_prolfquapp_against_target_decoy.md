# Review: protect prolfquapp against target+decoy FASTAs

Date: 2026-06-30
Reviewed plan: `TODO_protect_prolfquapp_against_target_decoy.md`

## Summary

I agree with the main direction: the root cause belongs upstream in annotation construction, and `ProteinAnnotation`
should guarantee one annotation row per protein ID before downstream report joins or `SummarizedExperiment` export see
the data. Moving the WU347806-specific FASTA dedup hotfix into a single `ProteinAnnotation` invariant is the right
shape.

The plan still needs tightening before implementation. The biggest risks are staging: Stage 1 removes the current REV
API and summary fields before all live consumers are enumerated, and the policy that "quant keeps decoys" conflicts
with current annotation-driven `get_subset(clean())` filtering. Stage 4 also undercounts row-annotation joins,
especially in QC report code.

This was a static review against the current working tree. I did not run tests.

## Findings

### 1. High: Stage 1 under-scopes the removal of REV/decoy consumers

The plan says Stage 1 should remove the `REV` column, `annotate_decoys()`, the `decoys=` branch of `clean()`/`nr_clean()`,
and REV-derived summary fields, then "Update `ProteinDataPrep`'s `get_summary()`/`remove_cont_decoy` consumers"
(`TODO_protect...`:140-144). That consumer list is incomplete, and `ProteinDataPrep$get_summary()` is not the current
method name.

Current consumers include at least:

- `R/R6_ProteinAnnotation.R:302-348`: `annotate_decoys()`, `percentOfFalsePositives`, `NrOfProteinsNoDecoys`.
- `R/R6_ProteinDataPrep.R:59-82`: `cont_decoy_summary()` and `remove_cont_decoy()` call `annotate_decoys()`,
  `nr_clean()`, and `clean(decoys=...)`.
- `inst/application/CMD_DEA_V2.R:220-223` and `R/cmd_helpers.R:458-461`: IBAQ subset still calls
  `xd$protein_annotation$clean(..., decoys = remove_decoys)`.
- `R/cmd_helpers.R:427-434`: logs `xd$protein_annotation$get_summary()`.
- `inst/application/DIANN/_Grp2Analysis_V2.Rmd:133-137`,
  `inst/application/CompoundDiscovery/Grp2Analysis_V2_R6.Rmd:142-146`, and
  `tests/testthat/Grp2Analysis_V2_R6.Rmd:147-151`: report prose reads the decoy summary fields and explains kept
  decoy sequences.
- `tests/testthat/test-R6_ProteinAnnotation.R:1-40`: current tests explicitly assert `REV`, `clean(decoys=...)`, and
  `nr_clean(decoys=...)` behavior.

Recommendation: split Stage 1 into two explicit substeps. First land the uniqueness/decoy-resolution invariant while
keeping compatibility shims or updating all direct callers. Then retire the report-facing REV fields and rewrite the
report text in the same commit that removes them. Otherwise Stage 1 is likely to produce a partially broken CLI/report
surface even if `ProteinAnnotation` itself is correct.

### 2. High: "LFQData keeps decoys" is not true with the current annotation-driven filtering path

The policy says the protein annotation is decoy-free while quant keeps decoys when no decoy pattern is given
(`TODO_protect...`:30-37). Current execution filters the quant data through the annotation:

- `R/R6_ProteinDataPrep.R:78-82`: `self$lfq_data_peptide <- self$lfq_data_peptide$get_subset(self$rowAnnot$clean(...))`.
- `inst/application/CMD_DEA_V2.R:220-223`: IBAQ uses `xd$lfqdata$get_subset(xd$protein_annotation$clean(...))`.
- `R/cmd_helpers.R:458-461`: helper path has the same IBAQ subset.

If `ProteinAnnotation` becomes decoy-free and `clean()` returns only annotation rows, any quant decoy with a distinct
protein ID will be dropped by this subset even if `remove_decoys` is retired. That contradicts the stated quant-side
policy and makes the accepted limitation at `TODO_protect...`:175-180 more consequential: distinct decoys are not only
undetected, they may be silently filtered depending on whether they remain in `row_annot`.

Recommendation: define the replacement semantics for `remove_cont_decoy()` and the IBAQ subset before removing
`clean(decoys=...)`. If quant decoys must be kept, the filter can only remove contaminants from quant rows, not "all
rows absent from a decoy-free annotation". If quant decoy retention is deferred with `prolfqua`, soften the current
policy text so it does not promise behavior this PR will not deliver.

### 3. High: Stage 4's join inventory is incomplete

The plan says there are about seven annotation-to-quant join sites and Stage 4 will "migrate the 7 sites"
(`TODO_protect...`:69-79,149-151). The current tree has more row-annotation joins than that, especially in QC:

- `R/R6_QC_Abundances.R:55-59`, `82-86`, `113-117`, `150-154`, `241-245`: left joins with
  `protein_annotation$row_annot`.
- `R/R6_QC_Abundances.R:194-199`, `226-230`, `467-471`: inner joins with `protein_annotation$row_annot`.
- `R/R6_QC_Abundances.R:122-126` and `250-257`: additional inner joins after annotation has already been attached,
  which can still drop rows if child-count data is incomplete.
- The plan correctly lists `R/R6_DEAnalyse.R:198-207`, `R/aggregation_IBAQ.R:83-87`, and
  `R/R6_DEAReportGenerator.R:230-266`.

Recommendation: replace the Stage 4 "7 sites" wording with an explicit inventory generated from the current source.
For each join, record the keeper side, explicit key, expected cardinality, and whether missing annotation is allowed.
Do not apply a generic helper to sample-level joins such as `resultList$annotation` in `R6_DEAReportGenerator`; the
helper should be scoped to row/protein annotation joins.

### 4. High: Stage 2's reader map is too coarse, and MaxQuant has duplicate implementations

The plan correctly identifies inconsistent `full_id` semantics, but the implementation map needs to be more precise:

- DIANN and BGS already pass `full_id = "fasta.id"` (`R/preprocess_DIANN.R:306-315`,
  `R/preprocess_BGS_default.R:197-205`).
- MaxQuant passes `full_id = "protein_Id"` in `R/preprocess_MaxQuant.R:438-446`. There is also a second
  `preprocess_MQ_peptide()` definition in `R/tidyMS_MaxQuant.R:367`, so both files need handling or one needs to be
  retired.
- FragPipe PSM passes `full_id = "protein_Id"` in `R/preprocess_FP_PSM.R:599-607`.
- MSstats is mixed: `preprocess_MSstats_FPDIA()` uses `full_id = "fasta.id"` at `R/preprocess_MSstats.R:173-181`, while
  the second MSstats path uses `full_id = "protein_Id"` at `R/preprocess_MSstats.R:280-288`.
- MaxQuant calls `get_annot_from_fasta(fasta_file)` without forwarding the configured `pattern_decoys`
  (`R/preprocess_MaxQuant.R:423`; same pattern in `R/tidyMS_MaxQuant.R:441`), so removing FASTA-level decoy filtering
  changes behavior there even before `full_id` is fixed.

Recommendation: make Stage 2 a table of reader functions, not just reader families. For each function, verify the join
key, the resulting `pID`, the raw prefixed column, and `cleaned_ids`. Update duplicated MaxQuant code in one deliberate
step.

### 5. Medium: `.is_decoy_id()` must explicitly ignore empty configured patterns

The incident driver is an empty `REVpattern`. `read_BF_yamlR6()` converts empty YAML patterns to `NULL`
(`R/R6_AppConfiguration.R:470-473`), but direct reader calls and current tests still pass `pattern_decoys = ""`
(`tests/testthat/test-get_annot_from_fasta-dedup.R:26`). If the new detector blindly unions an empty configured pattern
with defaults, `grepl("", ids)` matches every ID.

Recommendation: specify and test that `.is_decoy_id(ids, pattern = NULL)` and `.is_decoy_id(ids, pattern = "")` both use
only the built-in defaults, not an empty regex. Add an explicit test for this because it is the exact failure mode that
started this work.

### 6. Medium: `get_annot_from_fasta()` will not be a pure parser if raw-name dedup remains implicit

The plan says `get_annot_from_fasta()` should "return all rows (incl. decoys)" and let `ProteinAnnotation` own
decoy/duplicate resolution (`TODO_protect...`:118-123). The current function still removes duplicate raw FASTA record
names before parsing at `R/get_annot_from_FASTA.R:128`:

```r
fasta <- fasta[!(duplicated(names(fasta)))]
```

That may be acceptable for exact duplicate FASTA entries, but it is still a uniqueness decision outside
`ProteinAnnotation`.

Recommendation: decide explicitly. If this raw-name dedup stays, describe the parser contract as "all unique raw FASTA
records" and log the number removed. If the goal is truly that `ProteinAnnotation` owns all duplicate policy, move this
dedup or its accounting into the authority path too.

### 7. Medium: Stage 1 logging should distinguish "seen" from "dropped"

The plan wants construction logs for duplicated IDs, decoys dropped, and tiebreak-resolved counts (`TODO_protect...`:103-104).
Given the accepted within-duplicate-only limitation, a single "decoys dropped" count can hide important cases:

- decoy-like rows seen in duplicate groups and dropped,
- decoy-like rows seen outside duplicate groups and left untouched,
- duplicate groups resolved by `sp|`/`tr|`,
- duplicate groups resolved only by keep-first.

Recommendation: log at least those four counts. That makes the known limitation visible in production logs and avoids a
false sense that a target+decoy FASTA was fully cleaned when only duplicate groups were cleaned.

### 8. Medium: The temporary SE tryCatch needs an explicit final disposition

The working tree currently wraps `make_SummarizedExperiment()`/Quarto export in `tryCatch` in both
`inst/application/CMD_DEA_V2.R` and `inst/application/CMD_DEA_CD.R`. The plan says the WU347806 hotfix is a placeholder
to be folded in and retired, but Stage 6 only says "remove the placeholder hotfix" (`TODO_protect...`:153-155).

Recommendation: state whether the non-fatal SE export behavior is temporary or a permanent product decision. If the
root-cause invariant is mandatory, keeping SE export optional can mask regressions in the very path that exposed this
bug. If optional SE export is desired for operational reasons, add a targeted test that root-cause duplicate annotation
still fails earlier or is logged as a data invariant violation rather than only being swallowed at the export boundary.

### 9. Low: Status wording is stale against the current worktree

The plan says "Implementation not started" (`TODO_protect...`:3-5), but the current `prolfquapp` worktree already has
uncommitted changes in `R/get_annot_from_FASTA.R`, `R/R6_ProteinAnnotation.R`, `inst/application/CMD_DEA_V2.R`,
`inst/application/CMD_DEA_CD.R`, plus `TODO_fix_se_failure.md` and `tests/testthat/test-get_annot_from_fasta-dedup.R`.

Recommendation: clarify that the full streamline is not started, while a temporary WU347806 hotfix is already present
in the working tree.

## Test Additions I Would Require

- `ProteinAnnotation` duplicate group: target plus `REV_` decoy with same `pID`; forward row is kept and annotation is
  unique.
- Empty configured pattern: `pattern_decoys = ""` does not match all rows; defaults still catch `REV_`.
- Tiebreak duplicate group: `sp|ACC|...` beats `tr|ACC|...` when neither is a decoy.
- Fallback duplicate group: unrecognized duplicate rows still collapse deterministically and log keep-first.
- Distinct decoy limitation: a decoy with no forward twin is either explicitly retained with a warning/log count or the
  plan documents that it remains deferred.
- Reader-level fixtures for DIANN, BGS, MaxQuant, FP_PSM, and both MSstats paths proving `full_id` is the raw prefixed
  ID after preprocessing.
- Join-helper cardinality tests: duplicate annotation errors before join; missing annotation preserves quant/result row
  count where the plan says it should.
- SE regression: target+decoy FASTA no longer duplicates `diff_exp_analysis` and `make_SummarizedExperiment()` succeeds
  without relying on the SE tryCatch.

## Suggested Plan Edits Before Implementation

1. Amend Stage 1 to keep the uniqueness invariant separate from public/report API removal, or enumerate every caller and
   template that will be updated in the same commit.
2. Define the replacement for `remove_cont_decoy()` and IBAQ subset semantics so "annotation decoy-free" does not
   accidentally mean "quant decoys are dropped".
3. Replace the Stage 4 join count with a source-derived table of row-annotation joins.
4. Expand Stage 2 into a per-function reader matrix, including both MaxQuant definitions and both MSstats paths.
5. Add explicit empty-pattern behavior to the `.is_decoy_id()` design.
