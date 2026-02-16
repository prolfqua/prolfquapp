# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

Prolfquapp is an R package that provides a command-line interface (CLI)
for differential expression analysis (DEA) of proteomics and
metabolomics data. It wraps the `prolfqua` R package to streamline
analysis workflows in core facilities. The package processes outputs
from various quantification software (DIA-NN, MaxQuant, FragPipe, etc.)
and generates dynamic HTML reports with quality control plots and
statistical results.

## Common Development Commands

### R Package Development

``` bash
# Install package from source (development mode)
R CMD INSTALL .

# Build package
R CMD build .

# Check package (R CMD check)
R CMD check prolfquapp_*.tar.gz

# Run tests
Rscript -e "devtools::test()"
# or
Rscript -e "testthat::test_dir('tests/testthat')"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-FP_DDA_combined_protein.R')"

# Generate documentation (roxygen2)
Rscript -e "devtools::document()"

# Build vignettes
Rscript -e "devtools::build_vignettes()"
```

### Command-Line Tools (User-Facing)

The package installs 4 main CLI scripts in `inst/application/bin/`:

``` bash
# Copy shell scripts to working directory
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"

# 1. Generate annotation template from quantification outputs
./prolfqua_dataset.sh -i data_dir/ -s DIANN -d annotation.xlsx

# 2. Generate QC report
./prolfqua_qc.sh -i data_dir/ -p ProjectName -O ordername -w WorkunitName -d annotation.xlsx -s DIANN -o output_dir

# 3. Create YAML configuration file
./prolfqua_yaml.sh -y config.yaml

# 4. Run differential expression analysis
./prolfqua_dea.sh -i data_dir/ -d annotation.xlsx -y config.yaml -w AnalysisName -s DIANN

# All scripts support --help flag
./prolfqua_dea.sh --help
```

## Architecture

### R6 Classes (Object-Oriented Core)

The package uses R6 classes to manage state and orchestrate analysis:

- **ProlfquAppConfig** (`R6_AppConfiguration.R`): Central configuration
  manager
  - Contains nested R6 objects: ProcessingOptions, ProjectSpec,
    ExternalReader
  - Manages directory structure, normalization methods, FDR thresholds
  - Serializes to/from YAML for reproducible analyses
- **AnnotationProcessor** (`R6_AnnotationProcessor.R`): Experimental
  design handler
  - Reads annotation files (CSV/TSV/XLSX)
  - Validates required columns: file/raw, group/experiment,
    subject/bioreplicate, control
  - Automatically generates contrasts from control column (C/T
    designation)
  - Supports paired and unpaired designs
- **ProteinAnnotation** (`R6_ProteinAnnotation.R`): Protein metadata
  manager
  - Parses FASTA headers to extract gene names and annotations
  - Identifies contaminants (pattern: `^zz|^CON`) and decoys (pattern:
    `^REV`)
  - Provides filtering via `clean()` method
- **DEAnalyse** (`R6_DEAnalyse.R`): Complete DEA pipeline
  - Aggregates peptides to proteins (medpolish, lmrob, topN)
  - Transforms data (VSN, quantile, robscale, or none)
  - Builds 4 model types:
    - m1_linear: Protein-level moderated t-tests
    - m2_missing: Imputation-based analysis
    - m3_merged: Combined linear and missing results
    - m4_glm: Logistic models for missingness patterns
  - Computes contrasts with FDR correction
- **QC_generator** (`R6_QC_Abundances.R`): QC report generator
  - Computes IBAQ (intensity-based absolute quantification)
  - Generates protein abundance distributions
  - Renders R Markdown reports (QC_ProteinAbundances.Rmd, QCandSSE.Rmd)

### Software Preprocessing Pipeline

The preprocessing system uses a plugin architecture
(`preprocess_software.R`):

``` r
prolfqua_preprocess_functions <- list(
  "DIANN" = list(get_files, preprocess, extra_args, dataset),
  "MAXQUANT" = list(...),
  "FP_TMT" = list(...),
  # etc.
)
```

Each plugin handles: - File discovery (`get_files()`) - Data parsing and
filtering (`preprocess()`) - Software-specific parameters
(`extra_args`) - Annotation template generation (`dataset()`)

**Supported software:** - DIA-NN: `report.tsv` or `diann-output.tsv`
(q-value filtering) - MaxQuant: `proteinGroups.txt`, `peptides.txt`,
`evidence.txt` (LFQ intensities) - FragPipe: `psm.txt` (TMT, DIA, LFQ
workflows, PeptideProphet filtering) - MSstats: MSstats-formatted
input - Biognosys (BGS): BGS report format - MZMine: Metabolomics
feature tables

### Command-Line Scripts

Located in `inst/application/`:

- **CMD_MAKE_DATASET.R**: Generates annotation template from raw files
- **CMD_QUANT_QC.R**: Runs QC analysis and generates HTML reports
- **CMD_MAKE_YAML.R**: Creates YAML configuration files
- **CMD_DEA.R**: Main DEA pipeline (preprocess → aggregate → model →
  report)
- **CMD_DEA_Normalize.R**: Normalization-focused variant

All scripts use `optparse` for argument parsing and `logger` for
structured logging.

### Report Generation

R Markdown templates render dynamic HTML reports:

- \*\*\_Grp2Analysis_V2.Rmd\*\*: Main DEA report (methods, QC plots,
  volcano plots, results tables)
- \*\*\_Grp2Analysis_V2_metabo.Rmd\*\*: Metabolomics-specific variant
- \*\*\_DiffExpQC.Rmd\*\*: Differential expression QC
- **GenericQC/QC_ProteinAbundances.Rmd**: Protein abundance
  distributions
- **GenericQC/QCandSSE.Rmd**: Sample size estimation

Reports receive data via `params$GRP2` (list structure with results,
configurations, annotations).

### Output Formats

Multi-format outputs support different downstream tools:

- **HTML reports**: Human-readable results with plots
- **XLSX files**: Excel-compatible tables with multiple sheets
- **Parquet files**: Efficient storage of normalized data matrices
- **RDS files**: R-native serialization (SummarizedExperiment objects)
- **Rank files (.rnk)**: GSEA-compatible format
- **SummarizedExperiment**: Bioconductor-compatible for iSEE/OmicsViewer

## Key Dependencies

- **prolfqua**: Core statistical methods (LFQData class, modeling,
  contrasts, visualization)
- **arrow**: Parquet file I/O
- **SummarizedExperiment**: Bioconductor output format
- **optparse**: CLI argument parsing
- **logger**: Structured logging
- **rmarkdown/bookdown**: Dynamic report generation
- **seqinr**: FASTA parsing
- **vsn**: Variance-stabilizing normalization
- **UpSetR**: Missing data visualization

## Important Patterns

### Configuration-Driven Execution

All analysis parameters are stored in YAML files. R6 classes
serialize/deserialize configurations for reproducibility.

### Plugin System

New quantification software can be added by extending
`prolfqua_preprocess_functions` with:

``` r
"NEW_SOFTWARE" = list(
  get_files = function(data_dir) { ... },
  preprocess = function(xd) { ... },
  extra_args = list(...),
  dataset = function(xd) { ... }
)
```

### File-Based Pipeline

Sequential workflow with loose coupling: 1. Dataset creation →
annotation file 2. QC generation → HTML/XLSX reports 3. YAML creation →
config file 4. DEA execution → results directory

Each step reads files from previous steps, enabling modularity and
restart capability.

### Dynamic Function Resolution

The package uses string-encoded function names in configurations and
resolves them dynamically via
[`getFromNamespace()`](https://rdrr.io/r/utils/getFromNamespace.html),
allowing runtime customization of aggregation and transformation
methods.

## Testing

Tests are in `tests/testthat/` and follow naming convention `test-*.R`:

- `test-FP_DDA_combined_protein.R`: FragPipe DDA workflow
- `test-FP_DIA.R`: FragPipe DIA workflow
- `test-FP_TMT.R`: FragPipe TMT workflow
- `test-MQ_DDA.R`: MaxQuant DDA workflow

Tests use `testthat` framework and verify: - Output file generation
(.html, .xlsx, .rnk) - Directory structure creation - Expected number of
result files

## Core Facility Integration

The package integrates with B-fabric LIMS via standardized naming: -
`project_ID`: Project identifier - `order_ID`: Order identifier -
`workunit_ID`: Work unit identifier

Result directories are auto-organized by date and parameters (e.g.,
`DEA_2025-11-12_VSN_medpolish/`).
