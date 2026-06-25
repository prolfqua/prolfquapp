# TODO: Refactor the data-reader family (key normalization + join diagnostics + shared body)

**Status:** implementation plan. Approved scope so far: improve this TODO and start with the narrow DIANN
normalization / join-diagnostic work before any broad reader rewrite.
**Scope:** `prolfquapp/R/preprocess_*.R` + `R/tidyMS_MaxQuant.R` + `R/preprocess_software.R` +
`R/get_annot_from_FASTA.R`.
**Triggers:** two code-review items (`../../TODO/TODO_prolfqua_prolfquapp_code_review.md`):
1. *"Reader file-name normalization regex is inconsistent annot-side vs data-side (DIA-NN)"* — see Problem 1.
2. *"Reader `inner_join`s drop unmatched samples with only an info log"* — see Problem 5.

Both produce the **same symptom** (samples vanish at the annotation↔quant `inner_join`) at the **same join
sites**, but have different root causes — (1) is a wrong join *key*, (5) is the join silently tolerating
*genuine* key mismatches. Investigating (1) showed the incomplete annotation-side regex is copy-pasted across
**six readers**, sitting on top of a ~250-line copy-pasted reader body — so the narrow DIA-NN bug is one
symptom of a reader family that needs consolidating. See also `handover_automatic_reader_selection.md`.

## Why now

Commit `ce70dff` ("Fix regex patterns that incorrectly truncate file names") already fixed the **data-side**
DIA-NN regex (added `.d$`, `.mzML$`, and backslash normalization) but left the **annotation-side** copy
untouched, and `8cc364f` ("Read native DIA-NN 2.x output (parquet / Run column)") added DIA-NN 2.x. The two
sides have drifted, and the same drift-prone inline regex exists in five more readers. Any future per-reader
edit re-opens the same silent-sample-drop class of bug. Consolidating now removes the conditions for it.

## How readers are wired today (architecture)

- **Registry** (`R/preprocess_software.R:9-102`): `prolfqua_preprocess_functions` maps software keys
  (`DIANN`, `DIANN_PEPTIDE`, `FP_TMT`, `MAXQUANT`, `MSSTATS`, `MSSTATS_FP_DIA`, `BGS`, `MZMINE`, `SIM`,
  `DUMMY`, + `_PEPTIDE` variants) to a list of **string-encoded** function names and args:
  `get_files`, `preprocess`, `extra_args` (e.g. `"list(q_value = 0.01, hierarchy_depth = 1)"`), `dataset`.
  Resolved at runtime via `getFromNamespace(sub(...))` + `eval(parse(text = extra_args))`
  (`preprocess_software.R:139-151`).
- **Each `preprocess_*`** reads quant data, normalizes the sample key on the annotation side, joins
  annotation↔quant, builds `LFQData`, reads the FASTA via
  `get_annot_from_fasta()` (`R/get_annot_from_FASTA.R:101`, returns `fasta.id`, `fasta.header`,
  `proteinname`, `gene_name`, `protein_length`, `nr_tryptic_peptides`), builds a `ProteinAnnotation`, and
  returns `list(lfqdata, protein_annotation)`.
- **`get_*_files` / `dataset_template_*`** discover input files and build the annotation template.

## Problem 1 — Sample-key normalization is inconsistent, duplicated, and asymmetric

The sample key (the column the annotation↔quant `inner_join` matches on) is normalized **differently on the
two sides**, and the annotation-side regex is copy-pasted across six readers with the *incomplete* pattern
`"^x|\\.d\\.zip$|\\.raw$"` (missing `\\.d$`, `\\.mzML$`, and the backslash→slash step).

| Reader | Data-side key (verbatim source) | Annotation-side normalization | Match risk |
|---|---|---|---|
| DIANN | `Run`/`File.Name` → `gsub("^x\|\\.d\\.zip$\|\\.d$\|\\.raw$\|\\.mzML$","",basename(gsub("\\\\","/",x)))` (complete) | `gsub("^x\|\\.d\\.zip$\|\\.raw$","",basename(x))` (**incomplete**) | **`.d`/`.mzML`/Windows-path samples silently dropped** by `inner_join` (`preprocess_DIANN.R` join ~`:239`) |
| MaxQuant | column names `intensity.<name>` → strip `intensity.` (no basename, no ext-strip) | `tolower(gsub("^x\|\\.d\\.zip$\|\\.raw$","",basename(x)))` (**+ tolower**, reader-specific) | data side not lowercased here — matching depends on MQ column casing; `tolower` must be preserved/verified |
| MSstats / MSstats_FPDIA | `Run` **verbatim** | `gsub("^x\|\\.d\\.zip$\|\\.raw$","",basename(x))` | data key is whatever MSstats wrote; **adding `.d$`/`.mzML$` stripping could BREAK a join** if `Run` legitimately keeps that suffix |
| FP_PSM | TMT `channel` (column-name prefix stripped) | `gsub("^x\|\\.d\\.zip$\|\\.raw$","",basename(x))` → `raw.file` | join is `channel = file_name`; ext-strip rarely relevant to channel names |
| BGS | `R.FileName` **verbatim** | `gsub("^x\|\\.d\\.zip$\|\\.raw$","",basename(x))` | same verbatim-key risk as MSstats |
| MzMine | `datafile` (strip `datafile_`) | `basename(x)` only (**no ext-strip at all**) | yet another variant |

**Two distinct defects:**
1. **Drift (the named bug):** annotation side omits `.d$`, `.mzML$`, backslash-norm that the data side has.
   For DIANN this is a clear silent-drop bug and is **safe to fix** (data side already strips identically, so
   making the annotation side match is behavior-preserving except for the intended fix).
2. **Asymmetry is sometimes intentional:** MSstats/BGS take the data key **verbatim** (`Run`, `R.FileName`).
   There, the annotation normalization must reproduce *exactly what the software wrote* — so we **cannot
   assume** stripping `.d`/`.mzML` is correct without sample inputs. MaxQuant's `tolower` and MzMine's
   no-strip are likewise reader-specific. **This is why a blind one-regex-everywhere change is unsafe** and a
   per-reader, verified plan is needed.

## Problem 2 — ~250-line copy-pasted reader body (DRY)

DIANN, MaxQuant (`preprocess_MQ_peptide`), MSstats, MSstats_FPDIA, FP_PSM, BGS share an almost-identical
sequence: normalize annotation key → (filter `nr_peptides`) → set `config` columns (`file_name`,
`nr_children`, `ident_q_value`, `hierarchy[[protein_Id]]`, `hierarchy[[peptide_Id]]`, `set_response`,
`hierarchy_depth`) → `inner_join(annot, peptide)` → `prolfqua::setup_analysis` → `LFQData$new` →
`remove_small_intensities` → `get_annot_from_fasta` → `left_join(nrPEP, fasta_annot)` → `rename` →
`ProteinAnnotation$new` → `return list(lfqdata, protein_annotation)`. Only the **column mapping** differs:

| Reader | quant key | response | protein col | peptide col | fasta join key | cleaned_ids |
|---|---|---|---|---|---|---|
| DIANN | `raw.file` | `Peptide.Quantity` | `Protein.Group` | `Stripped.Sequence` | `Protein.Group.2 = proteinname` | `IDcolumn` |
| MaxQuant | `raw.file` | `peptide.intensity` | `leading.razor.protein` | `sequence` | `leading.razor.protein = fasta.id` | `proteinname` |
| MSstats | `Run`→file_name | `Intensity` | `ProteinName` | `PeptideSequence` | `ProteinName = fasta.id` | `proteinname` |
| MSstats_FPDIA | `Run`→file_name | `Intensity` | `ProteinName` | `PeptideSequence` | `ProteinName = proteinname` | `protein_Id` |
| FP_PSM | `channel` | `abundance` | `Protein` | `Peptide` | `Protein = fasta.id` | `proteinname` |
| BGS | `raw.file` | `FG.Quantity` | `PG.ProteinGroups` | `PEP.GroupingKey` | `Protein.Group.2 = proteinname` | `IDcolumn` |

## Problem 3 — Stringly-typed registry / dynamic dispatch

`extra_args` are stored as R source strings and `eval(parse())`-ed; functions resolved via string surgery on
`pkg::fun`. (Also flagged in the code review's "Dynamic dispatch" item.) Defeats static analysis, gives poor
errors. Lower priority than 1–2, listed for completeness.

## Problem 4 — Duplicated file-discovery helpers

`get_MQ_peptide_files` is defined **twice** (`preprocess_MaxQuant.R:309` and `tidyMS_MaxQuant.R:318`). The
six `get_*_files` share a discover-data + filter-`database*.fasta` + drop-`first-pass` pattern;
`dataset_template_*` each rebuild a `raw.file`/`Run`/`channel` unique list (some via the same readers, so
they inherit the data-side normalization — good — but should route through the same helper after refactor).

## Problem 5 — Silent partial-match drops at the annotation↔quant join

Every reader guards only `nr > 0` (at least one annotated file matched) and then `inner_join`s annotation
with quant data (`preprocess_DIANN.R` ~`:239`, `preprocess_MaxQuant.R:419`, `preprocess_MSstats.R:150,258`,
`preprocess_FP_PSM.R:578`, `preprocess_BGS_default.R:166`). If 8 of 10 annotated files match, the 2 unmatched
**silently vanish** and the analysis proceeds — surfaced only by an `info`-level log, easy to miss.

**Relationship to Problem 1.** Same symptom, same join sites, *different cause*: Problem 1 is a **wrong key**
(normalization drift makes keys that should match fail to); Problem 5 is the join **tolerating genuine key
mismatches** with no loud signal. They are complementary — fixing Problem 1 removes the
normalization-induced *false* drops, after which **any remaining drop is a real annotation↔data discrepancy**,
which is exactly what a Problem-5 diagnostic should surface. A correct key (1) makes the diagnostic (5)
meaningful instead of noisy.

**Intended-subsetting nuance (do NOT just `stop()`).** Per the maintainer note on the code-review item,
partial matching is *intentional*: the annotation is allowed to **subset** the data (exclude samples /
outliers, e.g. the annotation deliberately lists fewer runs than DIA-NN/Spectronaut produced). So the fix is
**not** "error on any partial match." It must distinguish:
- **intended subset** — annotation ⊂ data sample set → fine (this is how users drop samples); and
- **accidental drop** — annotated files that *should* match but don't (key mismatch, typo, format drift) →
  surface loudly (`log_warn`, or `stop()` only above a configurable tolerance).

Concretely: before/after the join, compare the distinct matched-sample count against the **annotation** count
(not the data count) and warn when annotated rows fail to match; report which annotated files were dropped.
Extra quant files that are not listed in the annotation are normal intended subsetting and should not warn by
themselves. The final home for this single, consistent check is the shared reader body, but it should be
introduced earlier as a small join/diagnostic helper so the safety fix is not blocked by the larger factory
refactor.

## Proposed refactor — staged, lowest-risk first

**Stage 1 — unify the DIANN key normalization and add the join diagnostic skeleton (resolves the urgent
code-review items without waiting for the large refactor).**
- Add one internal helper, e.g. `normalize_raw_file(x)`, implementing the **complete** data-side logic:
  cross-platform backslash handling (reuse/extend `normalize_path()` at `R/utils.R:75`) → `basename()` →
  `gsub("^x|\\.d\\.zip$|\\.d$|\\.raw$|\\.mzML$","",.)`. Home: `R/utils.R` (reusable) or `preprocess_DIANN.R`
  if kept DIANN-local. Internal (no `@export`) to keep API surface minimal — unless cross-package reuse
  (PTMreaders has its own copies) argues for exporting.
- **DIANN: apply now** by replacing both inline normalizers with the helper. This must preserve the current
  data-side behavior exactly; the intended behavior change is on the annotation side, where `.d`, `.mzML`, and
  Windows-style paths should now normalize the same way as quant keys.
- Add a small internal helper for the annotation-vs-quant sample check, e.g.
  `diagnose_sample_join(annotation_keys, quant_keys, matched_keys, context)`. In Stage 1 wire it into DIANN
  only, warning on annotated keys that did not match and allowing extra quant-only keys.
- **MaxQuant / MSstats / MSstats_FPDIA / FP_PSM / BGS / MzMine: per-reader verification gate.** For each,
  confirm the **data-side** key format before changing its annotation normalization:
  - readers whose data key is **verbatim** (`Run`, `R.FileName`) → must match exactly what the software
    writes; only adopt the shared helper if real sample keys are confirmed extension-free (else keep
    reader-specific behavior). Preserve MaxQuant's `tolower` and MzMine's no-strip unless proven safe to drop.
  - Decide the `^x` prefix quirk (strips a leading literal "x" from any name) — preserve as-is unless shown to
    misfire.

**Stage 2 — extract the shared body** into a factory, e.g.
`build_lfq_from_peptides(annot, quant, col_map, fasta_file, join_spec, patterns, hierarchy_depth, ...)`,
where `col_map`/`join_spec` carry the per-reader differences from the table above. Each `preprocess_*`
becomes a thin adapter: read+filter software data, build `col_map`, call the factory. (This is the code
review's `build_lfq_from_peptides` DRY item.) Move the Stage-1 join diagnostic helper into this factory so
all readers get the same annotation-missing warning without reimplementing it six times.

**Stage 3 (optional, lower priority)** — replace string-encoded `extra_args`/`getFromNamespace` with actual
lists + function references (or a typed registry).

**Stage 4** — collapse the duplicate `get_MQ_peptide_files`; parameterize `get_*_files` and route
`dataset_template_*` through the shared key helper.

## Implementation report — 2026-06-24

Stage 1 has been partially implemented for **DIANN only**.

Code changed:

- `R/preprocess_DIANN.R`
  - Added internal `.normalize_raw_file(x)`.
  - The helper strips the same wrappers on data-side and annotation-side keys:
    leading `x`, `.d.zip`, `.d`, `.raw`, `.mzML`, and Windows backslashes before `basename()`.
  - `diann_read_output()` now uses `.normalize_raw_file()` instead of the inline complete data-side regex.
    This is intended to preserve existing data-side behavior exactly.
  - `preprocess_DIANN()` now uses `.normalize_raw_file()` for annotation keys instead of the old incomplete
    annotation-side regex. This is the intentional behavior change: annotation values such as
    `C:\runs\sample.d` and `/mnt/data/control.mzML` now match DIA-NN quant keys `sample` and `control`.
  - Added internal `.diagnose_sample_join(annotation_keys, quant_keys, matched_keys, context)`.
  - `preprocess_DIANN()` captures annotation keys and quant keys before the join, then calls the diagnostic
    after `inner_join()`.
  - The diagnostic logs `log_warn` for annotation-only keys that did not match the joined data.
  - Quant-only keys are returned by the diagnostic but are not warned about, because those represent allowed
    user subsetting by annotation.

Tests changed:

- `tests/testthat/test-preprocess_DIANN-native2x.R`
  - Added a unit test for `.normalize_raw_file()` covering `.raw`, `.d`, `.mzML`, `.d.zip`, leading `x`,
    Windows backslashes, and plain names.
  - Added a diagnostic test proving annotation-only keys are reported and quant-only keys are allowed.
  - Added a diagnostic test proving quant-only keys do not produce log output.
  - Updated the existing `preprocess_DIANN()` native DIA-NN 2.x test so the annotation uses
    `C:\runs\sample.d` and `/mnt/data/control.mzML`, while the parquet `Run` values remain `sample` and
    `control`. This pins the original bug class directly.

Validation run:

- `air format R/preprocess_DIANN.R tests/testthat/test-preprocess_DIANN-native2x.R`
- `Rscript -e "devtools::test(filter = 'preprocess_DIANN-native2x|preprocess_DIANN-empty')"` → passed,
  `17` OK.
- `Rscript -e "devtools::test(filter = 'run_qc_preprocess')"` → passed, `9` OK, with two pre-existing SIM
  fixture warnings from `ProteinAnnotation$new()`.
- `git diff --check` → clean.

Validation not run:

- The installed-CLI integration fixture was not run after this code change. The integration harness uses the
  installed package script, so it requires refreshing the local install first. In this checkout, the integration
  `make install` target also installs the sibling `prolfqua` repo, which currently has unrelated dirty changes.
  Run the installed-CLI regression after deciding whether to install only `prolfquapp` or first stabilize the
  sibling repo state.

Remaining Stage 1 work:

- Run the small real installed-CLI regression on `integration_test/fixtures/diann_wu345302` after refreshing
  the installed package.
- Mine `/Users/wolski/projects/anndata_bridge` for broader real DIA-NN path/extension fixtures.
- Decide whether `.normalize_raw_file()` and `.diagnose_sample_join()` should stay DIANN-local or move to a
  shared helper file before applying them to additional readers.
- Do not apply the shared normalizer to MaxQuant, MSstats, MSstats_FPDIA, FP_PSM, BGS, or MzMine until their
  data-side key semantics have been verified from real fixtures.

## Constraints & risks

- `preprocess_*`, `diann_read_output`, `get_*_files`, `dataset_template_*` are **`@export`ed** and consumed by
  the CLI (`inst/application/CMD_*.R`), by `prolfquasaint` (deprecated redirect), and conceptually mirrored in
  `prolfquappPTMreaders` (which keeps its **own** copies of the regex — note `^x|.d.zip$|.raw$` there uses
  **unescaped** `.` = "any char"). Do not change public signatures; coordinate PTMreaders separately.
- The verbatim-key readers (MSstats/BGS) are the real trap: over-normalizing the annotation side can *create*
  silent drops. Stage-1 must be evidence-based per reader, not a blind sweep.
- Use `/Users/wolski/projects/anndata_bridge` as a broad real-data validation pool for reader edge cases
  (DIA-NN 2.x/native exports, `.d`/`.mzML` suffixes, path variants, and sample-name mismatch cases). Keep the
  smaller `integration_test/fixtures/diann_wu345302` fixture as the routine regression fixture.
- 120-char / 2-space lint applies to code; `make document` after any roxygen/signature change; NAMESPACE is
  generated.

## Verification (when implemented)

- Existing reader tests: `tests/testthat/test-preprocess_DIANN-empty.R`, `test-preprocess_DIANN-native2x.R`,
  `test-FP_DIA.R`, `test-run_qc_preprocess.R` must stay green.
- Add unit tests for `normalize_raw_file()` covering `.raw`/`.d`/`.mzML`/`.d.zip`/`^x`/Windows-backslash/plain.
- Add a per-reader **key-consistency** test: feed filenames with `.d`/`.mzML`/backslash and assert the
  annotation-side key equals the data-side key (so the `inner_join` retains all samples) — the durable
  regression guard against future drift.
- Synthetic `.d`/`.mzML` DIA-NN report + annotation → assert no samples dropped by the join.
- **Problem 5:** quant data listing files absent from the annotation (intended subset) → no warning;
  annotation listing files absent from the quant data (accidental drop) → `log_warn` fires (or `stop()` above
  the tolerance), and the reported dropped-file list is correct.
- Use `integration_test/fixtures/diann_wu345302` for the small real CLI regression. Because the integration
  harness runs the installed CLI, refresh the installed package first (`make install` / integration
  `make install`) before validating `CMD_DEA_V2.R` or reader changes through that path.
- Mine `/Users/wolski/projects/anndata_bridge` for broader real-reader fixtures once the narrow DIANN
  regression is green, especially cases involving native DIA-NN 2.x output, AnnData-derived exports,
  Windows-style paths, `.d`, and `.mzML`.

## Suggested order

1. Stage 1 for **DIANN only**: helper-based key normalization, DIANN join diagnostic, focused tests.
2. Validate with synthetic DIANN cases plus the small `diann_wu345302` installed-CLI regression.
3. Use `anndata_bridge` fixtures to choose the per-reader normalization policy for MaxQuant, MSstats,
   MSstats_FPDIA, FP_PSM, BGS, and MzMine before changing them.
4. Stage 2 (`build_lfq_from_peptides`) — extract the shared body and move the join diagnostic into that one
   place for all readers.
5. Then Stage 4 (file-discovery dedupe), then Stage 3 (typed registry).
