# FGCZ Report Homogenisation — Quarto POC

## Context

FGCZ produces HTML reports from multiple omics pipelines (RNA-seq via ezRun/SUSHI, proteomics via prolfquapp). Currently these have inconsistent styling, structure, and no machine-readable metadata. The goal is a unified, branded Quarto report format that:

1. Works for any omics type (RNA-seq, proteomics, spatial, etc.)
2. Has a recognisable FGCZ identity (banner + clean theme, polished later)
3. Includes agent-readable metadata as a sidecar JSON file
4. Follows the homogenisation spec (self-contained, tabbed, callout boxes, lightboxes, downloadable PDFs of all figures in `plots/`)

Both test datasets are `SummarizedExperiment` objects — this is the shared container.

## Test data

- **RNA-seq**: `/srv/gstore/projects/p38719/o38990_DESeq2_E185--over--E165_2026-03-04--12-05-19/E185--over--E165/` (output.rds, deResult.rds, param.rds) — two-group DESeq2, has GO/ORA/GSEA in metadata
- **Proteomics**: `https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342063/3106962.rds` (downloaded to `/tmp/3106962.rds`) — 6877 proteins, 12 samples, 3 groups (Manual/OT/SandwichOnly), two contrasts

## Design decisions

- **All tables use `DT::datatable()`** — no `kable()`/`kableExtra`. Consistent interactive tables throughout.
- **Agent metadata is a sidecar file** — `report-metadata.json` written alongside the HTML, not embedded in it.
- **Data loading uses `qs2::qs_read()` with `readRDS()` fallback** — via a helper that auto-detects format.
- **Homogenised file names** — both reports follow the pattern `diffExp_<omics>.qmd` (e.g. `diffExp_rnaseq.qmd`, `diffExp_proteomics.qmd`).

## Deliverables

All files go in `inst/templates/quarto/`.

### 1. `_fgcz-report.yml` — shared Quarto metadata
Defines format defaults both reports include:
- `format: html` with `embed-resources: true`
- `theme: cosmo` (minimal for now)
- FGCZ header include
- `lightbox: true`
- `code-fold: true`, `code-tools: true`
- `fig-dpi: 300`
- Cross-reference config
- Execute defaults (echo: false, warning: false, message: false)

### 2. `fgcz_header_quarto.html` — Quarto-compatible header
Slimmed version of current header:
- FGCZ banner (inline base64)
- Enrichr JS function
- External-link targeting JS

### 3. `_helpers.R` — shared R functions
- `fgczLoadData(name)` — loads `.qs2` or `.rds` files (auto-detect)
- `fgczSavePlot(plot, name, width, height)` — saves to `plots/` as PNG + PDF with informative names
- `fgczAgentMetadata(...)` — writes `report-metadata.json` sidecar file
- `fgczDatatable(df, ...)` — thin wrapper around `DT::datatable()` with FGCZ defaults (compact, striped, scrollX, filter)

### 4. `diffExp_rnaseq.qmd` — RNA-seq differential expression report
Quarto conversion of twoGroups.Rmd. Structure:

```
::: {.panel-tabset}

# Settings
  - Callout box with params (DT::datatable)
  - Timestamp

# Results Summary
  - Reference, feature counts
  - Significant counts (DT::datatable)
  - Result xlsx download link
  - Live report link

# Significant Genes {.panel-tabset}
  ## Between-group comparison
    - Scatter plot (ggplot2 static, saved to plots/)
    - Volcano (p-value, saved to plots/)
  ## Advanced
    - Volcano (FDR)
    - P-value histogram
  ## Intra-group
    - Per-group scatter plots

# Clustering
  - Callout with thresholds
  - Heatmap (pheatmap)
  - GO cluster table link

# Enrichment {.panel-tabset}
  ## Enrichr
  ## ORA {.panel-tabset}
    - Overview + per-ontology tabs with DT::datatable
  ## GSEA {.panel-tabset}
    - Overview + per-ontology tabs with DT::datatable

# Technical Bias
  - Fisher test (DT::datatable)

# Input Dataset
  - DT::datatable

# Session Info

:::
```

Key changes from Rmd:
- `{.tabset}` -> `:::{.panel-tabset}` fenced divs
- **All kable() calls replaced with DT::datatable()**
- All plots get `fgczSavePlot()` to `plots/` with PDF variant
- Callout boxes for settings/thresholds/method descriptions
- Lightbox on figure chunks
- Agent metadata written as sidecar JSON at end
- Data loading via `fgczLoadData()` (qs2 with rds fallback)

### 5. `diffExp_proteomics.qmd` — Proteomics DE report
Parallel structure, adapted for proteomics SE:

```
::: {.panel-tabset}

# Settings
  - Callout box: project, order, formula, contrasts (DT::datatable)
  - B-Fabric links from metadata$bfabric_urls

# Results Summary
  - Protein counts, significant counts per contrast (DT::datatable)
  - Result xlsx download

# Quality Control {.panel-tabset}
  ## Missing Values
    - NA heatmap (present/absent)
  ## Abundance Distribution
    - Density plots (raw vs normalised)
  ## CV Analysis
    - Per-group CV violin
  ## PCA
    - PCA from transformedData assay

# Differential Expression {.panel-tabset}
  ## Overview
    - Callout: method, formula, contrasts (DT::datatable)
    - Significant counts
  ## Volcano
    - Per-contrast volcano plots (saved to plots/)
  ## Results Table
    - DT::datatable with all proteins

# Significant Proteins
  - Filtered DT::datatable (FDR < threshold)
  - Heatmap of significant proteins
  - UpSet across contrasts

# Technical QC
  - Model diagnostics (imputed vs moderated counts, DT::datatable)

# Input Dataset
  - Sample annotation (DT::datatable)

# Session Info

:::
```

### 6. `report-metadata.json` — agent/robot sidecar file
Written by `fgczAgentMetadata()` at render time:

```json
{
  "@context": "https://fgcz.ch/report-schema/v1",
  "reportType": "differential-expression",
  "omicsType": "rnaseq|proteomics",
  "organism": "...",
  "comparison": [
    { "name": "X_vs_Y", "sample": "X", "reference": "Y", "method": "DESeq2|prolfqua" }
  ],
  "thresholds": { "pValue": 0.05, "log2FC": 0.5, "fdr": 0.05 },
  "counts": { "totalTested": N, "significantUp": N, "significantDown": N },
  "topHits": [ { "id": "...", "name": "...", "log2FC": N, "fdr": N } ],
  "files": { "resultTable": "result.xlsx", "plots": "plots/", "report": "diffExp_rnaseq.html" },
  "qcFlags": [],
  "generatedAt": "ISO-8601",
  "softwareVersions": { "R": "...", "quarto": "..." }
}
```

## File structure

```
inst/templates/quarto/
  _fgcz-report.yml            # shared Quarto metadata defaults
  fgcz_header_quarto.html     # banner + JS
  _helpers.R                   # shared R functions
  diffExp_rnaseq.qmd          # RNA-seq DE report
  diffExp_proteomics.qmd      # Proteomics DE report
```

## Implementation order

1. **Shared infrastructure**: `_fgcz-report.yml`, `fgcz_header_quarto.html`, `_helpers.R`
2. **RNA-seq report**: `diffExp_rnaseq.qmd` — convert from twoGroups.Rmd, test with p38719 data
3. **Proteomics report**: `diffExp_proteomics.qmd` — build from p23078 SE, parallel structure
4. **Agent metadata**: implement `fgczAgentMetadata()`, wire into both reports
5. **Verify**: render both, check self-containment, visual consistency

## Critical files to reference

- Current template: `/misc/GT/analysis/peter/ezRun/inst/templates/twoGroups.Rmd`
- Current header: `/misc/GT/analysis/peter/ezRun/inst/templates/fgcz_header.html`
- Current CSS: `/misc/GT/analysis/peter/ezRun/inst/templates/fgcz.css`
- Rendering pipeline: `/misc/GT/analysis/peter/ezRun/R/util.R` (`makeRmdReport`, `ezLoadRobj`)
- Report helpers: `/misc/GT/analysis/peter/ezRun/R/reports.R` (`makeResultFile`, `makeCountResultSummary`)
- RNA-seq test data: `/srv/gstore/projects/p38719/o38990_DESeq2_E185--over--E165_2026-03-04--12-05-19/E185--over--E165/`
- Proteomics test data: `/tmp/3106962.rds`
- Quarto system: `/opt/quarto/bin/quarto` v1.8.25

## Verification

1. Copy test data to a working directory under `inst/templates/quarto/`
2. `quarto render diffExp_rnaseq.qmd` — produces self-contained HTML
3. `quarto render diffExp_proteomics.qmd` — produces self-contained HTML
4. Both have FGCZ banner, same tabbed layout, callout boxes, DT tables
5. `plots/` directory created with PDF + PNG for every figure
6. `report-metadata.json` present alongside each HTML with valid JSON
7. Lightboxes work on figures
8. Reports open correctly without internet (self-contained)
