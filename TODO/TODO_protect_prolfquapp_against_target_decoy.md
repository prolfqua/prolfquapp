# Protect prolfquapp against target+decoy FASTAs — design & implementation plan

**Date:** 2026-06-30
**Status:** **IMPLEMENTED** (working tree; not committed). Stages 1-5 + 7 done; full testthat suite
**239 pass / 0 fail**; validated on the real WU347806 FASTAs. Stage 6 (reader builder
consolidation) **deferred** — see note at the end. The temporary WU347806 hotfix has been folded
in / removed: the `get_annot_from_fasta` `sp|` dedup and the keep-first guard are replaced by the
real design; the non-fatal SE `tryCatch` is removed (SE export is mandatory again).
**Companions:** `prolfqua/TODO/TODO_revpattern_handling.md` (deferred core work **+ the full
REV-API removal/reintroduction inventory**); `TODO_fix_se_failure.md` (original crash root-cause).
**Reviewed:** decisions below incorporate `REVIEW_protect_prolfquapp_against_target_decoy.md`.

## Implementation summary (what landed)

- **Stage 1** — `R/R6_ProteinAnnotation.R`: added internal `.detect_decoy_ids()` (defaults ∪
  configured pattern; `""`/`NULL`/`"a^"` → defaults only) and `.resolve_unique_protein_ids()`
  (within-duplicate: decoy-drop → `sp|` tiebreak → keep-first; logs 4 counts); `initialize()`
  calls the resolver; added `get_rev_pattern()`; removed `annotate_decoys()` + the `decoys=` arms of
  `clean()`/`nr_clean()` + the REV summary fields; `clean()` now does **pattern-gated** decoy
  removal. Cascade updated: `R6_ProteinDataPrep.R` (`cont_decoy_summary`, `remove_cont_decoy`),
  `cmd_helpers.R` + `CMD_DEA_V2.R` IBAQ subsets, and the 4 source report templates. Tests rewritten
  in `test-R6_ProteinAnnotation.R`.
- **Stage 2** — no change: verified every reader's `full_id` already references a prefixed-id column
  (by-`fasta.id` readers → `protein_Id` is prefixed; by-`proteinname` readers → `fasta.id`). The
  earlier "fix to fasta.id" would have broken MaxQuant/FragPipe/MSstats.
- **Stage 3** — `get_annot_from_FASTA.R` is now a pure parser (removed the raw-name dedup, the decoy
  pre-filter, and the `sp|` collapse). `test-get_annot_from_fasta-dedup.R` rewritten.
- **Stage 4** — added `.join_annotation()` (right-join, explicit `by`, asserts unique annotation);
  migrated the 2 `get_annotated_contrasts` + 5 `prep_result_list` joins. QC/IBAQ joins left as-is
  (already `left_join`/explicit-`by`, safe under the unique-annotation invariant). Tests in
  `test-report_helpers.R`.
- **Stage 5** — deleted `R/tidyMS_MaxQuant.R` (stale near-duplicate of `preprocess_MaxQuant.R`);
  migrated the unique `dataset_template_MAXQUANT` into `preprocess_MaxQuant.R`. NAMESPACE unchanged.
- **Stage 6** — **deferred** (stretch): a single shared reader builder would refactor 5+ readers
  with thin/no test coverage (no MaxQuant test exists) for a DRY-only gain; the consolidation
  primitive `build_protein_annot()` already exists for future use.
- **Stage 7** — SE export reverted to mandatory (loud & fatal) in `CMD_DEA_V2.R`/`CMD_DEA_CD.R`.
  Validated on the real WU347806 FASTAs (parser passes 58,333 records; resolver collapses every
  target+decoy collision to the unique forward). **Not yet run:** the full 538 MB end-to-end rerun
  + production image rebuild.

## What actually happened (one paragraph)

The reference was a standard FGCZ **target+decoy** UniProt DB
(`p24227_canislupus_..._d_...fasta`): every forward entry (`tr|`, `sp|`, contaminant `zz|`) has a
`REV_`-prefixed decoy twin → 54,251 records, 27,124 forward + 27,124 decoy + 3 specials. The run
was configured with `REVpattern: ''` (empty), so `get_annot_from_fasta()` never stripped the
decoys. The accession extractor `gsub(".+\\|(.+)\\|.*", "\\1", fasta.id)` discards everything
before the first `|` — including the `REV_` prefix — so `tr|ACC|…` and `REV_tr|ACC|…` both collapse
to `ACC`. Result: ~50% duplicated protein IDs (`mean(duplicated) = 0.4999724`), which fanned out
through the `multiple = "all"` joins and finally aborted `make_SummarizedExperiment()` at
`column_to_rownames()` ("duplicate 'row.names' are not allowed"). `protein.fas` was a strict subset
of the species FASTA and contributed nothing — it was not a cause. (Note: the raw-name dedup at
`get_annot_from_FASTA.R:128` could not catch this — it dedups full raw ids, which were unique; the
collision is at the extracted-accession level, one step later.)

## Why this didn't surface years ago

1. **Decoys are normally stripped early** — with a correct `REVpattern`, `get_annot_from_fasta()`
   removes them before IDs can collapse. This run had an empty pattern.
2. **`ProteinAnnotation` never enforced ID uniqueness** — nothing checked the invariant.
3. **Old outputs tolerated duplicates silently** (XLSX/HTML don't set row names → just bloat). The
   crash only appeared once the **new** SummarizedExperiment export routed those IDs through
   `column_to_rownames()`, which requires unique names.

## Decisions locked in (Q&A + review resolution)

1. **Full streamline** — one shared path, implemented in stages.
2. **`prolfqua` core untouched now**; rev-pattern→`LFQData` handoff deferred (its own TODO).
3. **`ProteinAnnotation` is the single authority for ID uniqueness.** Resolve duplicates **within
   duplicate-ID groups only** (for now): decoy-drop → `sp|`>`tr|` tiebreak → **keep-first
   fallback** (always yields one row per ID, even when `full_id` has no usable prefix).
4. **Retire the REV machinery atomically** — remove the `REV` column, `annotate_decoys()`,
   the `decoys=` branch of `clean()`/`nr_clean()`, and the REV-derived summary fields, updating
   every caller / report template / test **in one coordinated change**. Full inventory lives in
   `prolfqua/TODO/TODO_revpattern_handling.md`.
5. **`get_rev_pattern()` getter** exposes the detected pattern (future handoff to `LFQData`).
6. **Quant-side decoy removal is pattern-gated**, not annotation-driven: a configured
   `pattern_decoys` → decoys removed from quant before DEA; **no** configured pattern → quant keeps
   them (decoy-proportion QC deferred). This replaces the dual-purpose `clean(decoys=)`.
7. **SE export stays mandatory** — failures are loud & fatal; the `tryCatch` hotfix is removed once
   the invariant lands.
8. **Contaminants stay** (`CON` / `annotate_contaminants` / `clean(contaminants=)` /
   `remove_cont`); the `remove_decoys` boolean is retired (`remConDec` still drives contaminants).
9. **Single uniqueness authority** — also remove the FASTA parser's raw-name dedup
   (`get_annot_from_FASTA.R:128`); the parser returns every record verbatim.

## Decoy policy (reconciled)

- **The annotation carries one row per protein ID.** Decoys that **collide with a forward**
  (the duplicate case) are dropped during de-duplication. Decoys that appear as their **own
  distinct ID** are *not* removed for now (within-duplicate-only — see Known limitation).
- **Quant decoy removal is gated on the configured pattern** (decision 6).
- **No annotation→quant join may drop or multiply quant/result rows.**

## Current code reality (reader map + review)

- **Decoy handling is split across three layers** today: `get_annot_from_fasta()` pre-filters
  decoys *only if* `pattern_decoys` is non-empty (the WU347806 gap) and collapses duplicate
  accessions (hotfix); `ProteinAnnotation` flags `REV` lazily in `annotate_decoys()` and
  keep-first-dedups. The refactor makes `ProteinAnnotation` the sole authority.
- **`full_id` (the column decoy patterns match) is inconsistent** — see reader matrix.
- **The `remConDec` YAML param drives both** contaminant and decoy removal; only the contaminant
  half survives. Empty `REVpattern` is already coerced to `NULL` in `read_BF_yamlR6()`
  (`R6_AppConfiguration.R:471`), but readers/tests still pass `pattern_decoys = ""`
  (e.g. `preprocess_MSstats.R:94`), so the detector must treat `""`/`NULL` identically.
- **Annotation↔quant join sites (source-derived, regenerate at Stage 4):**
  - `ProteinAnnotation$initialize` — `left_join(distinct(pID), row_annot)` then the
    `nr_children_experiment` **`inner_join`** (can drop pIDs).
  - `ProteinDataPrep$remove_cont_decoy` — `get_subset(clean(...))` (R6_ProteinDataPrep.R:78-82).
  - `DEAnalyse$get_annotated_contrasts` — `inner_join(rowAnnot$row_annot, …, multiple="all")` ×2
    (R6_DEAnalyse.R:198,203).
  - `aggregation_IBAQ` — `inner_join` (line 83).
  - `R6_QC_Abundances` — left joins (55-59, 82-86, 113-117, 150-154, 241-245), inner joins
    (122-126, 194-199, 226-230, 250-257, 467-471).
  - `DEAReportGenerator$prep_result_list` — 5× `inner_join(row_annot, …, multiple="all")` with
    **implicit `by=`** (230-266).
  - `cmd_helpers.R` / `CMD_DEA_V2.R` — IBAQ `get_subset(clean(decoys=…))`.

### Reader matrix (`full_id` + decoy forwarding)

| reader (fn / file) | `full_id` today | nrPep↔fasta join key | action |
|---|---|---|---|
| DIANN (`preprocess_DIANN.R`) | `fasta.id` ✅ | `proteinname` | none |
| BGS (`preprocess_BGS_default.R`) | `fasta.id` ✅ | — | verify |
| FP_DIA (`preprocess_MSstats_FPDIA`) | `fasta.id` ✅ | `proteinname` | note default `pattern_decoys=""` |
| MSstats (`preprocess_MSstats`) | `protein_Id` ⚠️ | `fasta.id` | **fix → `fasta.id`** |
| MaxQuant (`preprocess_MaxQuant.R`) | `protein_Id` ⚠️ | `fasta.id` | **fix**; also forwards no `pattern_decoys` to `get_annot_from_fasta` (L423) |
| FragPipe (`preprocess_FP_PSM.R`) | `protein_Id` ⚠️ | `fasta.id` | **fix → `fasta.id`** |
| MzMine / CompoundDisc | feature id | — | no-op (`pattern_decoys=NULL`, no decoys) |

**Duplicate MaxQuant impl:** `preprocess_MQ_peptide` is defined twice — `preprocess_MaxQuant.R:349`
(canonical) and `tidyMS_MaxQuant.R:367` (legacy, shadows by load order). `tidyMS_MaxQuant.R` is a
near-complete stale duplicate; **only `dataset_template_MAXQUANT` (L345) is unique to it** and is
wired into the registry. Plan: migrate `dataset_template_MAXQUANT` into `preprocess_MaxQuant.R`,
diff the duplicated bodies to confirm canonical, then delete `tidyMS_MaxQuant.R`.

## Target design

### 1. Shared decoy-id detector (prolfquapp, internal)

`.is_decoy_id(ids, pattern = NULL)` → logical. Built-in **anchored** default prefixes (`^REV_`,
`^rev_`, `^DECOY`, `^decoy_`, `^XXX_`, `^reverse_`, `^##`) **unioned** with a configured pattern.
**`pattern = NULL` or `pattern = ""` use the defaults only** — never `grepl("", …)` (the exact
incident driver). Returns the matched pattern so `ProteinAnnotation` can expose it.

### 2. `ProteinAnnotation` = single uniqueness + decoy authority

After building `row_annot`, within each duplicated-`pID` group: (1) drop rows where
`.is_decoy_id(row[[full_id]], pattern_decoys)` (keep forward); (2) tiebreak `sp|`>`tr|` on
`full_id` for UniProt ids; (3) **keep-first fallback**. Always one row per `pID`. Store + expose
via `get_rev_pattern()`. Change the `nr_children_experiment` join to a quant-preserving `left_join`.
**Log four counts:** duplicated IDs seen, decoys dropped (within dups), `sp|`-resolved,
keep-first-resolved. Remove `REV` column, `annotate_decoys()`, `decoys=` branches, REV summary
fields; keep all contaminant machinery.

### 3. Standardize `full_id` = the prefixed identifier

`full_id` must be the raw, prefixed id (`fasta.id`). Fix MaxQuant / FragPipe / MSstats (regular
path) per the matrix; verify their quant `protein_Id` is or isn't prefixed against real inputs.

### 4. `get_annot_from_fasta()` → pure parser

Remove the decoy pre-filter, the `sp|` collapse (hotfix), **and the raw-name dedup (line 128)**.
Returns every record as read; `ProteinAnnotation` owns all uniqueness.

### 5. One quant-preserving join helper (row/protein-annotation joins only)

`join_annotation(quant_or_result, annotation, by = pID)` → `left_join` keeping the quant/result
side, **explicit `by`**, asserting unique annotation (cannot multiply). Migrate the row/protein
joins above. **Do not** apply it to sample-level joins (e.g. `resultList$annotation`).

### 6. Quant-side decoy removal (replaces `clean(decoys=)`)

`ProteinDataPrep$remove_cont_decoy()` and the IBAQ subsets: remove contaminants always
(`remove_cont`), and remove decoys **only when a `pattern_decoys` is configured**, using
`get_rev_pattern()` against the quant ids. No pattern → quant keeps decoys.

### 7. (Stretch) reader builder consolidation

Extract `compute_nr_peptides_per_protein()` + `build_protein_annotation()` so readers stop
duplicating count→join→rename→`$new`.

## Implementation plan (staged; each stage = commit + tests)

- **Stage 1 — `ProteinAnnotation` authority (atomic REV retirement).** `.is_decoy_id()` (with
  empty-pattern handling); within-duplicate resolution + 4-count logging; `get_rev_pattern()`;
  `nr_children` join → `left_join`; remove REV column/`annotate_decoys()`/`decoys=`/REV summary
  fields **and** update every caller, the report templates, and the tests in the same change (use
  the inventory in `prolfqua/TODO/TODO_revpattern_handling.md`). Define quant-side decoy removal
  (§6) so "annotation decoy-free" never silently drops quant decoys.
- **Stage 2 — standardize `full_id`** (reader matrix); per-reader tests decoys are detected/dropped.
- **Stage 3 — `get_annot_from_fasta()` pure parser** (drop pre-filter, `sp|` collapse, line-128
  dedup); update its tests; verify `prolfquappPTMreaders`' `fasta_annot_early` path still holds.
- **Stage 4 — quant-preserving join helper**; migrate the row/protein joins; explicit `by` in
  `prep_result_list`; tests that no quant/result row is dropped or multiplied; verify report code
  tolerates `NA`-annotation rows.
- **Stage 5 — drop legacy MaxQuant** (`tidyMS_MaxQuant.R`): migrate `dataset_template_MAXQUANT`,
  diff bodies, delete.
- **Stage 6 (stretch) — reader builder consolidation** (§7).
- **Stage 7 — validate & retire hotfix.** End-to-end rerun of recovered `WU347806input.zip`;
  confirm `SummarizedExperiment.rds` + `outputs.yml`; **remove the SE `tryCatch`** so a target+decoy
  FASTA succeeds with no safety net; update `TODO_fix_se_failure.md`.

## Tests to add (from review)

- Duplicate group: target + `REV_` twin same `pID` → forward kept, annotation unique.
- `pattern_decoys = ""` does not match all rows; defaults still catch `REV_`.
- Tiebreak: `sp|ACC|…` beats `tr|ACC|…` when neither is a decoy.
- Fallback: unrecognized duplicates collapse deterministically and log keep-first.
- Per-reader fixtures (DIANN, BGS, MaxQuant, FP_PSM, both MSstats paths) prove `full_id` is the raw
  prefixed id after preprocessing.
- Join-helper cardinality: duplicate annotation errors before join; missing annotation preserves
  quant/result row count.
- SE regression: target+decoy FASTA no longer duplicates `diff_exp_analysis`;
  `make_SummarizedExperiment()` succeeds **with the `tryCatch` removed**.

## Risks / watch

- **Inner→left on result joins** may surface `NA`-annotation rows previously dropped — but the
  primary safety is the unique-annotation invariant (Stage 1), which kills multiplication
  everywhere; the join switch is robustness. Verify report code tolerates `NA` annotation.
- **`cont_decoy_summary()`** (the real method name) + the report prose lose decoy fields — rewrite
  the four source report templates (vignette, test template, DIANN, CompoundDiscovery); leave the
  generated `prolfquasaint/inst/.../Inputs_WU_*` artifacts alone.
- **`full_id` fix** must be verified against real MaxQuant/FragPipe/MSstats outputs.

## Known limitation (accepted for now)

Within-duplicate-only detection does **not** catch a decoy that appears as its own distinct quant
ID (no forward twin) — acceptable for now (DIANN reports forward proteins only). Such a distinct
decoy will **remain** in both the annotation and the quant (since quant removal is pattern-gated and
the decoy isn't a duplicate). Eager whole-annotation detection + the decoy-proportion QC are the
future step in `prolfqua/TODO/TODO_revpattern_handling.md`. Logging makes this visible rather than
implying a target+decoy FASTA was fully cleaned.
