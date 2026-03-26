# Add DEA Results to AnnData

## Motivation

The preprocess methods can already create AnnData objects.
Our R6_DEAReportGenerator writes the results as Excel files, rank files and many more.

The idea: start adding the results to the AnnData object primarily.
Then have a second module function which generates from the AnnData file all the files we currently generate directly, including the HTML files and reports.

**Why this is a good move:**

1. **Single source of truth** — DEA results (fold changes, p-values, contrasts) live in one structured object rather than being scattered across Excel/RNK/RDS writers
2. **Decouples computation from presentation** — the analysis pipeline produces an AnnData; a separate export module reads it and generates HTML, Excel, rank files, etc. Both halves are independently testable
3. **Interoperability** — AnnData (h5ad) is readable from Python, opening the door to downstream tools (scanpy ecosystem, custom dashboards, `anndata_omics_bridge`)
4. **Reproducibility** — a single `.h5ad` file captures normalized data + results + metadata, easier to archive and re-render reports from


## The Pipeline Vision

```
                         CMD_DEA_V2 (extended)
                        /                     \
CMD_PREPROCESS    or  existing flags           export module
      |               (-i -d -s)                    |
  raw data + annot        |                    reads DEA .h5ad
  + fasta + software      |                         |
      |                   v                         v
      v             aggregation + DEA          XLSX, ORA, GSEA,
  peptide.h5ad ---→ + writes protein-level     HTML reports, SE,
  (-a flag)         dea_protein.h5ad           boxplots, etc.
```

1. **CMD_PREPROCESS** (new) — standalone CLI step: takes raw data + annotation + FASTA + software flag, writes peptide-level `.h5ad`
2. **CMD_DEA_V2** (extended) — accepts either `.h5ad` via new `-a` flag OR existing `-i`/`-d`/`-s` flags. Runs aggregation + DEA, writes protein-level `.h5ad` enriched with DEA results
3. **Export module** — reads the protein-level `.h5ad` and generates all output files


## Existing Infrastructure

We already have substantial AnnData support via **anndataR** (pure R, no reticulate):

| File | What it does |
|------|-------------|
| `R/preprocess_DIANN_anndata.R` | `preprocess_DIANN_anndata()` — preprocesses DIA-NN data and returns AnnData. `preprocess_anndata_from_lfq()` — generic LFQData+ProteinAnnotation to AnnData converter |
| `R/lfq_to_ann_data.R` | `anndata_from_LFQData()` — simpler LFQData-to-AnnData conversion (layers, obs, var) |
| `R/anndata_to_lfqdata.R` | `LFQData_from_anndata()` — **round-trip back** from AnnData to LFQData+ProteinAnnotation. Uses `uns$prolfquapp` namespace for lossless reconstruction. This is the entry point for prolfqua_dea to consume an .h5ad |

The `uns` namespace already has a schema (`schema_version: "1.0.0"`) with:
- `uns$prolfquapp` — analysis_configuration, protein_annotation, layer_names (for round-trip)
- `uns$exploreDE` — column role metadata (anndata_omics_bridge convention)
- `uns$X_layer_name` — name of primary intensity column


## Column Name Indirection (AnalysisConfiguration Pattern)

Nothing is hardcoded. All column and layer names are software-dependent and resolved through `AnalysisConfiguration`, which is serialized into `uns$prolfquapp$analysis_configuration`. The key fields:

| Config field | What it names | Example values |
|-------------|---------------|----------------|
| `workIntensity` | Intensity layer(s) — a stack, last = current response | `"peptide.intensity"`, `"log2_intensity"` |
| `nr_children` | Nr peptides/precursors per protein per sample (layer) | `"nr_peptides"`, `"nr_precursors"`, `"nr_psm"` |
| `ident_qValue` | Identification q-value (layer) | `"qValue"` |
| `ident_Score` | Identification score (layer, if present) | |
| `hierarchy` | Feature identity columns (var) | `list(protein_Id = "protein_Id")` |
| `factors` | Experimental design columns (obs) | `list(group = "group", subject = "subject")` |
| `fileName`, `sampleName` | Sample identity columns (obs) | |

`config$value_vars()` returns all value columns that become layers. `config$annotation_vars()` returns all sample annotation columns for obs. Consumers look up the indirection — never assume fixed names.

This is already implemented in `preprocess_anndata_from_lfq()` and `LFQData_from_anndata()`.


## Protein Annotation in Peptide-Level AnnData

**Problem:** At peptide level, `var` has one row per peptide, but `ProteinAnnotation` has one row per protein (far fewer rows). It doesn't fit as `var` (wrong cardinality) or `obs` (that's samples).

**Current solution** (in `preprocess_anndata_from_lfq()`, line 173): **left join** protein annotation onto peptide-level `var` by `pID`. Protein annotation columns (description, cleaned_ids, exp_nr_children, contaminant/decoy flags) are **repeated** for all peptides of the same protein.

This is pragmatic and works well because:
- It keeps `var` self-contained — every peptide row carries its protein context
- No extra slot needed
- Standard AnnData tooling can filter/group by protein columns directly
- The join key (`pID`) is stored in `uns$prolfquapp$protein_annotation$pID` so consumers know which column links peptides to proteins

The `ProteinAnnotation` field names are also stored in `uns$prolfquapp$protein_annotation` (full_id, description, cleaned_ids, exp_nr_children, pattern_contaminants, pattern_decoys) for round-trip reconstruction.


## Data Model: Protein-Level DEA AnnData

```
adata (samples x proteins)
|
+-- X                     # primary intensity matrix (name from config$get_response())
+-- layers                # ALL value columns from config$value_vars(), names are NOT fixed:
|   |                     #   - workIntensity stack (e.g. "peptide.intensity", "log2_intensity")
|   |                     #     the raw vs transformed distinction comes from the stack order
|   |                     #   - nr_children (e.g. "nr_peptides", "nr_precursors", "nr_psm")
|   |                     #   - ident_qValue (e.g. "qValue")
|   |                     #   - ident_Score, opt_rt, opt_mz (if present)
|   |                     # All names resolved via uns$prolfquapp$analysis_configuration
|   +-- <workIntensity[1]>  # e.g. raw intensity
|   +-- <workIntensity[2]>  # e.g. transformed intensity (= X)
|   +-- <nr_children>       # e.g. "nr_peptides" — per sample per protein
|   +-- <ident_qValue>      # e.g. "qValue" — per sample per protein
|   +-- ...
|
+-- obs                   # sample annotation (config$annotation_vars():
|                         #   fileName, sampleName, factor columns, normValue)
|
+-- var                   # protein annotation (from ProteinAnnotation$row_annot)
|                         #   hierarchy columns (config$hierarchy_keys())
|                         #   description, cleaned_ids, full_id
|                         #   <exp_nr_children> (experiment-wide summary, column name from
|                         #     uns$prolfquapp$protein_annotation$exp_nr_children)
|
+-- varm                  # one entry per contrast
|   +-- "condition_A_vs_B"   # data.frame: diff, statistic, p.value, FDR (rows = proteins)
|   +-- "condition_C_vs_D"   # data.frame: diff, statistic, p.value, FDR (rows = proteins)
|   +-- ...
|
+-- uns
    +-- prolfquapp
    |   +-- schema_version    # "2.0.0" (bumped — now includes DEA results)
    |   +-- source_software   # "DIANN", "MAXQUANT", etc.
    |   +-- analysis_configuration  # full AnalysisConfiguration fields (for round-trip)
    |   +-- protein_annotation      # ProteinAnnotation field names (for round-trip)
    |   +-- layer_names             # from config$value_vars()
    |   +-- dea                     # NEW: DEA-specific metadata
    |       +-- formula             # model formula string
    |       +-- default_model       # which model was used (e.g. "m1_linear")
    |       +-- FDR_threshold       # FDR cutoff used
    |       +-- diff_threshold      # fold-change cutoff used
    |       +-- summary             # model summary table
    |       +-- contrast_names      # list of contrast names (keys into varm)
    |       +-- contrasts           # contrast definitions (name -> formula string)
    +-- exploreDE
        +-- column_roles    # anndata_omics_bridge convention (already present)
```


## Data Model: Peptide-Level Preprocessing AnnData (already implemented)

```
adata (samples x peptides)
|
+-- X                     # primary intensity (config$get_response())
+-- layers                # all config$value_vars() — same indirection as protein-level
|
+-- obs                   # sample annotation (same as protein-level)
|
+-- var                   # peptide annotation + protein annotation LEFT-JOINED
|                         #   hierarchy columns at full depth (protein_Id, peptide_Id, ...)
|                         #   protein annotation duplicated across peptides of same protein
|                         #   (description, cleaned_ids, exp_nr_children, contaminant/decoy flags)
|                         #   join key: uns$prolfquapp$protein_annotation$pID
|
+-- uns
    +-- prolfquapp        # same structure as protein-level (schema 1.0.0)
    +-- exploreDE         # column roles
```


## What DEAReportGenerator Currently Produces

From `R6_DEAReportGenerator.R`, the `write_DEA_all()` method generates:

| Output | Source in code | Maps to AnnData |
|--------|---------------|-----------------|
| XLSX (14 sheets: annotation, raw/norm abundances, contrasts, stats, etc.) | `prep_result_list()` + `writexl::write_xlsx()` | Derivable from obs, var, X, layers, varm, uns |
| ORA gene lists (.txt per contrast direction) | `.write_ORA()` from `annotated_contrasts_signif` | Derivable from varm (filter by FDR + diff sign) |
| GSEA rank files (.rnk per contrast) | `.write_GSEA()` from `annotated_contrasts` | Derivable from varm (extract statistic column) |
| HTML DEA report | `render_DEA()` with Grp2Analysis_V2_R6.Rmd | Needs the AnnData + Rmd template |
| HTML QC report | `render_DEA()` with DiffExpQC_R6.Rmd | Needs the AnnData + Rmd template |
| Per-protein boxplot PDFs | `write_protein_boxplots()` | Needs abundance data from X/layers |
| SummarizedExperiment (.rds) | `make_SummarizedExperiment()` | Derivable from the same AnnData |

All tabular outputs (XLSX, ORA, GSEA, SE) are derivable from the AnnData structure above.


## Implementation Plan

### Phase 0: CMD_PREPROCESS + extend CMD_DEA_V2 to accept AnnData

**CMD_DEA_V2.R and all existing infrastructure continue to work exactly as-is.** We add a new preprocessing CLI step and extend CMD_DEA_V2 to also accept `.h5ad` as an alternative input.

1. **New CMD_PREPROCESS.R** — standalone CLI script that takes:
   - `-i data_dir/` (quantification output)
   - `-d annotation.xlsx` (sample annotation)
   - `-f protein.fasta` (FASTA file)
   - `-s SOFTWARE` (DIANN, MAXQUANT, FP_TMT, etc.)
   - `-o output.h5ad`

   Calls the appropriate `preprocess_*()` → `preprocess_anndata_from_lfq()` → `anndataR::write_h5ad()`.
   The generic converter already exists; only DIA-NN has a dedicated `preprocess_DIANN_anndata()` but the pattern generalizes trivially since all preprocessors return `list(lfqdata, protein_annotation)`.

2. **New shell wrapper `prolfqua_preprocess.sh`** (alongside existing `prolfqua_dea.sh` etc.)

3. **Extend CMD_DEA_V2.R** with a new `-a input.h5ad` flag as an alternative to the existing `-i`/`-d`/`-s` flags:
   - When `.h5ad` is provided, uses `LFQData_from_anndata()` (already exists) to reconstruct LFQData + ProteinAnnotation
   - When raw data flags are provided, existing preprocessing path runs as before
   - Rest of the DEA pipeline is identical in both cases

After this phase there are two input paths into CMD_DEA_V2, both working:
- **Existing:** `CMD_DEA_V2 -i data_dir/ -d annotation.xlsx -s DIANN ...` (preprocessing + DEA in one step, as today)
- **New:** `CMD_PREPROCESS → peptide.h5ad → CMD_DEA_V2 -a peptide.h5ad → results`

### Phase 1: DEA writes results back into AnnData

1. **New function `dea_results_to_anndata(deanalyse, source_software)`** — takes a completed `DEAnalyse` object and returns a protein-level AnnData:
   - `X` = primary intensity from `config$get_response()` (last entry in workIntensity stack)
   - `layers` = all `config$value_vars()` pivoted wide — **no hardcoded names**, everything driven by AnalysisConfiguration
   - `obs` = sample annotation from `config$annotation_vars()`
   - `var` = protein annotation from `deanalyse$rowAnnot$row_annot` + hierarchy columns from `config$hierarchy_keys()`
   - `uns$prolfquapp$analysis_configuration` = full AnalysisConfiguration state — the indirection layer that tells consumers which layer is `nr_children`, which is the response, etc.
   - `varm` = one entry per contrast from `deanalyse$contrast_results[[default_model]]$get_contrasts()`, pivoted wide by contrast name
   - `uns$prolfquapp$dea` = formula, thresholds, summary, contrast definitions
2. **Add `write_anndata()` method to DEAReportGenerator** — calls `dea_results_to_anndata()` and writes `.h5ad` via `anndataR::write_h5ad()`
3. **Call from `write_DEA_all()`** — alongside existing writers, so nothing breaks
4. **Bump schema to 2.0.0** in `uns$prolfquapp$schema_version`

After this phase the full pipeline is: `preprocess → peptide.h5ad → CMD_DEA_V2 → dea_protein.h5ad`

### Phase 2: Export module that reads AnnData and produces all outputs

1. **New function `export_from_anndata(h5ad_path, outdir)`** — reads the `.h5ad` and generates:
   - XLSX with the same 14 sheets as `prep_result_list()`
   - ORA text files (filter varm by FDR/diff threshold)
   - GSEA rank files (extract statistic from varm)
   - SummarizedExperiment .rds
2. **Report rendering from AnnData** — either:
   - Convert AnnData back to DEAnalyse (via round-trip) and render with existing Rmd templates, or
   - Write new Rmd templates that accept AnnData directly (cleaner long-term, bigger effort)
3. **Deprecate direct writers** in DEAReportGenerator once the AnnData export path is validated

### Phase 3 (future): Python consumers

- `anndata_omics_bridge` / `anndata_proteomics_bridge` can read the `.h5ad` natively
- Build Python-side export tools, dashboards, or downstream analysis that consume the same file
- The `exploreDE` column roles namespace enables generic tooling without prolfquapp-specific knowledge
