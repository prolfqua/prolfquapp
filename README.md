[![DOI](https://zenodo.org/badge/642819797.svg)](https://doi.org/10.5281/zenodo.15845909)

# prolfquapp (ˌproʊˈlɛf.kə.ˌæp): Generating Dynamic DEA Reports using a command line interface to the prolfqua R Package

Read on JPR <https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00911> "prolfquapp ─ A User-Friendly Command-Line Tool Simplifying Differential Expression Analysis in Quantitative Proteomics"

*Prolfquapp* is a command-line interface to the [prolfqua](https://github.com/fgcz/prolfqua) R package ([doi](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00441)) for protein differential expression analysis. It preprocesses outputs from DIA-NN, MaxQuant, FragPipe, and Spectronaut, and generates HTML reports, Excel tables, rank files, and SummarizedExperiment objects for downstream tools such as [ExploreDE](https://github.com/fgcz/exploreDE).

![prolfquapp](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/Prolfqapp_Highlight.png?raw=true)

# Differential Expression Analysis Workflow with prolfquapp

After running your quantification software (DIA-NN, MaxQuant, FragPipe-TMT, FragPipe-DIA, or FragPipe-LFQ), place the quantification outputs and the `.fasta` file used for the search into a single `data_dir`.

Copy the shell scripts into your working directory:

``` bash
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

or using the [Docker](https://www.docker.com/products/docker-desktop/) container:

``` bash
prolfquapp_docker.sh R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

This places five scripts into your working directory:

``` bash
[1] "/<working_directory>/prolfqua_dea.sh"
[2] "/<working_directory>/prolfqua_yaml.sh"
[3] "/<working_directory>/prolfqua_qc.sh"
[4] "/<working_directory>/prolfqua_dataset.sh"
[5] "/<working_directory>/prolfqua_contrasts.sh"
```

On Linux, make them executable:

``` bash
chmod a+x prolfqua_*
```

All scripts support `--help`. All commands can be prefixed with `./prolfquapp_docker.sh` to run in the Docker container instead of a local R installation.

## Workflow Overview

1.  [Create Dataset](#1-create-dataset)
2.  [Generate Quality Control (QC)](#2-generate-quality-control-qc)
3.  [Generate prolfquapp YAML](#3-generate-prolfquapp-yaml)
4.  [Generate Contrast Definitions (optional)](#4-generate-contrast-definitions-optional)
5.  [Run Differential Expression Analysis](#5-run-differential-expression-analysis)

## 1. Create Dataset

Generate an experiment annotation template from the quantification output files.

-   **Input**: directory containing identification/quantification software outputs
-   **Output**: annotation file (CSV, TSV, or XLSX)

``` bash
./prolfqua_dataset.sh -i data_dir/ -s DIANN -d annotation.xlsx
```

The generated `annotation.xlsx` contains five columns:

-   `Relative.Path` / `Path` / `raw.file` / `channel` — file identifier (must be unique)
-   `name` — label used in tables and figures (must be unique)
-   `group` / `experiment` — main factor
-   `subject` / `bioreplicate` — blocking factor (optional; delete column if experiment is unpaired)
-   `control` — reference condition marker (`C` = control, `T` = treatment) (optional)

The `raw.file` column is pre-filled from the input directory. Fill in the remaining columns before proceeding.

## 2. Generate Quality Control (QC)

Generate a QC report consisting of two HTML documents and an XLSX file.

-   **Input**: annotation file from step 1 and quantification output directory
-   **Output**: subfolder starting with `QC_` containing QC report and visualizations

``` bash
./prolfqua_qc.sh -i data_dir/ -p ProjectName -O ordername -w WorkunitName -d annotation.xlsx -s DIANN -o where_to_write_results
```

## 3. Generate prolfquapp YAML

Create a YAML configuration file with the DEA parameters.

-   **Output**: YAML configuration file

``` bash
./prolfqua_yaml.sh -y config.yaml
```

Edit the generated YAML file to set any additional parameters not exposed via the command line.

## 4. Generate Contrast Definitions (optional)

Add contrast information to the annotation file.

``` bash
# Single factor: adds CONTROL column (C = reference, T = rest)
./prolfqua_contrasts.sh annotation.xlsx --control WT -o annotation_with_control.xlsx

# Two factors: adds ContrastName/Contrast columns
./prolfqua_contrasts.sh annotation.xlsx --f1 treatment --f2 time -o annotation_with_contrasts.xlsx
```

## 5. Run Differential Expression Analysis

Run the DEA using the annotation and configuration files from the previous steps.

-   **Input**: quantification output directory, annotation file (step 1 or 4), YAML config (step 3)
-   **Output**: subfolder starting with `DEA_` containing HTML reports, Excel tables, rank files, and `SummarizedExperiment.rds`

``` bash
./prolfqua_dea.sh -i data_dir/ -d annotation.xlsx -y config.yaml -w NameOfAnalysis -s DIANN
```

## How to install

**Linux**

``` bash
export R_LIBS_SITE="/scratch/PROLFQUA/r-site-library/"
R --vanilla << EOF
.libPaths()
install.packages(c("remotes","seqinr", "prozor", "logger", "arrow"), repos = "https://stat.ethz.ch/CRAN/")
remotes::install_gitlab("wolski/prolfquadata", host="gitlab.bfabric.org")
remotes::install_github("fgcz/prolfqua", build_vignettes = TRUE, dependencies = TRUE)
remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
EOF
```

**Docker**

Download [prolfquapp_docker.sh](https://raw.githubusercontent.com/prolfqua/prolfquapp/refs/heads/master/inst/application/bin/prolfquapp_docker.sh) and use it as a prefix to any command (see above).

## ASMS poster: Streamlining Protein Differential Expression Analysis in Core Facilities

![prolfquapp_ASMS_poster](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/prolfquapp_PosterPNG.png?raw=true)

# Related software

-   Einprot <https://github.com/fmicompbio/einprot>
-   LFQAnalyst <https://analyst-suite.monash-proteomics.cloud.edu.au/apps/lfq-analyst/> and <https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst>
-   POMAShiny <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009148>
-   MSDap <https://github.com/ftwkoopmans/msdap>
