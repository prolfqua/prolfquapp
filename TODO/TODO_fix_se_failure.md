# Fix: WU347806 DEA failure — SummarizedExperiment export aborts on duplicate protein IDs

**Date:** 2026-06-30
**Status:** Root cause **FIXED at the source** by the target+decoy streamline (see
`TODO_protect_prolfquapp_against_target_decoy.md`). `ProteinAnnotation` now guarantees unique
protein IDs (decoy-aware within-duplicate resolution), so the `column_to_rownames` abort cannot
recur. The temporary non-fatal SE `tryCatch` hotfix has been **removed** — SE export is mandatory
again (loud & fatal). Validated: full testthat suite 239 pass / 0 fail; and on the real WU347806
FASTAs the new parser passes all 58,333 records through and the resolver collapses every
target+decoy collision to the unique forward (0 decoys kept). **Still pending:** the full 538 MB
end-to-end rerun (`make install` + run the recovered inputs) to confirm `SummarizedExperiment.rds`
+ `outputs.yml` on disk, and the production image rebuild.
**Workunit:** 347806 (`DEA_firth_C1_allTissuesVSTumor`), order/project 32824.

## Symptom

The DEA run produced all of the main deliverables (`DE_WU347806.xlsx`, `DE_WU347806.html`,
`QC_WU347806.html`, `lfqdata_normalized.parquet`, ORA/GSEA files, `IBAQ_.xlsx`) and then
**aborted** while writing the SummarizedExperiment. The prolfqua log
(`prolfqua_202606251108.log`, recovered from `WU347806input.zip`) stops at:

```text
INFO [2026-06-25 11:14:18] Writing summarized experiment.
```

with no further lines. Because the process died (uncaught R error → non-zero exit) before
`outputs.yml` was written, app-runner registered no output resource and B-Fabric marked the
workunit failed. The real R error went to stderr, not this log.

The abort point in the source is the **mandatory** SE/Quarto block in
[inst/application/CMD_DEA_V2.R:240-251](../inst/application/CMD_DEA_V2.R#L240-L251):

```r
logger::log_info("Writing summarized experiment.")
se_file <- file.path(reporter$resultdir, "SummarizedExperiment.rds")
SE <- reporter$make_SummarizedExperiment()   # <-- throws here
saveRDS(SE, file = se_file)
```

## Root cause (CONFIRMED)

A duplicated `protein_Id` in the protein annotation propagates into the per-contrast result
tables; `make_SummarizedExperiment()` then tries to set those `protein_Id`s as data.frame row
names, which base R forbids when they are not unique.

### The chain

1. **Two overlapping FASTA files were supplied** for this run:
   `p24227_canislupus_UP000002254_9615_d_20230511.fasta` (54,251 records) and
   `protein.fas` (4,082 `sp|` records). The species FASTA already lists **each accession
   roughly twice**, and `protein.fas` overlaps it.

2. [`get_annot_from_fasta()`](../R/get_annot_from_FASTA.R#L101) dedups only by the raw record
   name ([get_annot_from_FASTA.R:128](../R/get_annot_from_FASTA.R#L128)
   `fasta <- fasta[!(duplicated(names(fasta)))]`), then extracts the UniProt accession into
   `proteinname` ([:174](../R/get_annot_from_FASTA.R#L174)). It **detects** duplicate
   `proteinname` ([:178-189](../R/get_annot_from_FASTA.R#L178-L189)) but, because
   `REVpattern` was empty (`config.yaml` → `7|REVpattern: ''`), takes the `warning()` branch
   and **returns the annotation with duplicates intact**. (If a decoy pattern had been set it
   would have `stop()`-ed instead — a different, conflated guard.)

3. [`ProteinAnnotation$initialize`](../R/R6_ProteinAnnotation.R#L254-L262) builds `row_annot`
   as `distinct(protein_Id)` from the data and then **left-joins the FASTA annotation**:

   ```r
   self$row_annot <- dplyr::distinct(dplyr::select(lfqdata$data_long(), self$pID))
   if (!is.null(row_annot)) {
     self$row_annot <- dplyr::left_join(self$row_annot, row_annot, by = self$pID)
   }
   ```

   Because `row_annot` has duplicate `protein_Id`, this left-join **multiplies** every protein
   row. No `distinct()` follows, so `self$row_annot` now carries duplicate `protein_Id`.

4. [`prep_result_list()`](../R/R6_DEAReportGenerator.R#L230-L251) joins `row_annot` into the
   contrasts and matrices with `multiple = "all"`:

   ```r
   ctr <- dplyr::inner_join(ra$row_annot, contr_obj$get_contrasts(), multiple = "all")
   ```

   so `diff_exp_analysis` (and the wide matrices) inherit the duplicated `protein_Id` rows.
   This is also why **`DE_WU347806.xlsx` is a bloated 44 MB** — the table is duplicated.

5. [`make_SummarizedExperiment()`](../R/R6_DEAReportGenerator.R#L581-L656) splits
   `diff_exp_analysis` per contrast and calls
   [`column_to_rownames(..., var = "protein_Id")`](../R/R6_DEAReportGenerator.R#L637-L644):

   ```r
   for (i in names(diffbyContrast)) {
     row.data <- prolfquapp::column_to_rownames(diffbyContrast[[i]], var = rowname)
     ...
   }
   ```

   `column_to_rownames` ([report_helpers.R:10-16](../R/report_helpers.R#L10-L16)) does
   `rownames(res) <- <protein_Id vector>`. With duplicate values base R throws:

   ```text
   Error: duplicate 'row.names' are not allowed
   ```

   → uncaught → process aborts exactly at "Writing summarized experiment."

### Why earlier steps succeeded

The XLSX/HTML/Parquet writers never set row names, so duplicated rows are written without
error (just bloat). The matrices `mat.raw`/`mat.trans` come from `lfq_data$data_wide()` (which
has unique `protein_Id`), so the SE assays themselves are fine. The failure is specific to the
**rowData** assembly, which routes the duplicated `diff_exp_analysis` through
`column_to_rownames`.

### Why other workunits did not hit this

They were run with a single, non-duplicated FASTA, so `row_annot` had unique `protein_Id`.

## Reproduction (done — cheap, no full pipeline rerun needed)

Both ends of the chain reproduced locally against the real inputs from `WU347806input.zip`
(extracted into `inst/application/DIANN/`):

```r
# (1) the two real FASTAs -> duplicate accessions at the production rate
fa <- prolfquapp:::get_annot_from_fasta(c(species_fasta, "protein.fas"), pattern_decoys = "")
nrow(fa)                       # 54251
dplyr::n_distinct(fa$proteinname)   # 27127
mean(duplicated(fa$proteinname))    # 0.4999724   <-- matches log: 0.499972350740079

# (2) column_to_rownames aborts on a duplicated key (tibble input, as diff_exp_analysis is)
df <- tibble::tibble(protein_Id = c("P1","P1","P2"), diff = c(1,2,3))
prolfquapp::column_to_rownames(df, var = "protein_Id")
#> Error: duplicate 'row.names' are not allowed
```

The `0.4999724` duplicate fraction is byte-for-byte the value in the production log, so this is
the same condition. A full end-to-end rerun (538 MB `diann-output.tsv`, ~6 min) is optional
confirmation only.

## Fix

Per the repo policy (fix the root cause upstream, no bandaids), the annotation must yield
**one row per protein ID**. Two complementary changes:

### Primary (root cause) — dedup `proteinname` in `get_annot_from_fasta` — APPLIED

[R/get_annot_from_FASTA.R:178-205](../R/get_annot_from_FASTA.R#L178-L205). The old block
conflated "duplicate IDs" with "wrong decoy pattern" and otherwise returned duplicates with a
misleading warning. Now: when a decoy pattern **is** set, duplicates still `stop()` (forward and
reverse most likely share an accession — the pattern didn't match). When no decoy pattern is set,
the duplicates are real overlapping inputs, so we collapse to one row per `proteinname`,
**preferring the reviewed `sp|` entry** (per your call) by ordering `sp|` rows first
(`order(db != "sp")`, where `db <- sub("\\|.*$", "", fasta.id)`) before `!duplicated()`. Both
the warn (count + rate) and the collapse (unique count) are logged.

### Complementary (invariant guard) — enforce unique `protein_Id` in `ProteinAnnotation` — APPLIED

[R/R6_ProteinAnnotation.R:281-296](../R/R6_ProteinAnnotation.R#L281-L296). At the exit of the
constructor, after both annotation joins, any duplicated `protein_Id` rows are collapsed
(keeping the first) with a `log_warn`. This holds the "one row per protein ID" invariant
regardless of annotation source (proteinGroups, anndata `obs`, user tables), not just FASTA.

### Make SE/Quarto export non-fatal — APPLIED (per your go-ahead, with loud logging)

The live paths — [CMD_DEA_V2.R:240-272](../inst/application/CMD_DEA_V2.R#L240-L272) and
[CMD_DEA_CD.R:311-340](../inst/application/CMD_DEA_CD.R#L311-L340) — now wrap the
`make_SummarizedExperiment()` / `saveRDS()` / `render_quarto_se_report()` chain in `tryCatch`.
On failure it logs at **ERROR** level with the `conditionMessage` **and the failing call**
(`conditionCall`), then continues to `write_index_html`, archiving, zipping, and the `outputs.yml`
write — so a supplementary-export failure can no longer sink an otherwise-complete DEA, and the
cause stays in the run log.

Note: `cmd_helpers.R::write_dea_run_outputs()` (an extracted, already-`tryCatch`'d copy of the
output-writing pipeline) is currently **unused** — the CMD scripts duplicate the logic inline.
Consolidating the two onto that single helper is a worthwhile follow-up but was left out of scope
to keep this change minimal.

## Regression test — ADDED

[tests/testthat/test-get_annot_from_fasta-dedup.R](../tests/testthat/test-get_annot_from_fasta-dedup.R):

1. Two overlapping FASTAs (`tr|ACC1` + `sp|ACC1`, same accession) → asserts unique `proteinname`
   and that the **reviewed `sp|` copy is the one kept** (via its distinct `GN=` gene name).
2. A decoy pattern set + duplicate IDs still `stop()`s with `wrong decoys pattern` (guard preserved).

Full suite green via `devtools::load_all()`: **216 pass / 0 fail** (10 pre-existing benign
warnings, 15 CRAN-only skips).

## Validation — remaining (not yet run)

```bash
cd prolfquapp
make install        # reinstall as 2.2.8 (installed runtime is currently 2.2.7)
```

End-to-end: rerun the recovered inputs (`inst/application/DIANN/WU347806input.zip`) through
`CMD_DEA_V2.R` with `model=firth`, `config.yaml`, `dataset.csv`; confirm `SummarizedExperiment.rds`
is written and `outputs.yml` is produced. (If `quarto` is absent locally, the Quarto render now
logs an error and continues — `SummarizedExperiment.rds` + the rest still complete.) Then
rebuild/publish the prolfquapp image and rerun 347806, or salvage the already-generated result
files from the recovered zip.

## Open questions

- **Keep-first vs. prefer reviewed (`sp|`)** when collapsing duplicate accessions — confirm the
  preferred policy.
- **Why is the species FASTA doubled?** `p24227_canislupus_…fasta` already contains each
  accession ~twice. Worth flagging to whoever builds these DBs; prolfquapp should be robust
  regardless, but a doubled reference DB may indicate an upstream DB-prep bug.
- Do you also want the **non-fatal SE export** change applied now, or kept separate?
