# prolfquapp: Generating Dynamic DEA Reports using a command line interface to the prolfqua R Package

![prolfquapp](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/Prolfqapp_Highlight.png?raw=true)

Uses preprocessing and statistical models implemented in the R package [prolfqua](https://github.com/fgcz/prolfqua)
[doi.org/10.1021/acs.jproteome.2c00441](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00441)

## How to install

On Linux

```bash
export R_LIBS_SITE="/scratch/PROLFQUA/r-site-library/"
R --vanilla << EOF
.libPaths()
install.packages(c("remotes","seqinr", "prozor", "logger"), repos = "https://stat.ethz.ch/CRAN/")
remotes::install_gitlab("wolski/prolfquadata", host="gitlab.bfabric.org")
remotes::install_github("fgcz/prolfqua", build_vignettes = TRUE, dependencies = TRUE)
remotes::install_github("prolfqua/prolfquapp", dependencies = TRUE)
EOF
```


## How To use prolfquapp

prolfquapp is a command line tool. To use it open you shell (linux, mac), or command window (windows).
Change into the directory with the identification/quantification results coming from FragPipe, MaxQuant, DIA-NN, Spectronaut etc.
In the directory execute:

```bash
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

This will place the following four shell script files (linux), or bat files (windows) into your working directory:


```
[1] "/<working_directory>/prolfqua_dea.sh"
[3] "/<working_directory>/prolfqua_yaml.sh"
[4] "/<working_directory>/prolfqua_qc.sh"
[5] "/<working_directory>/prolfqua_dataset.sh"
```

On Linux give the executables LINUX permissions:

```
chmod a+x prolfqua_*
```

All scripts can be run with the option `--help`.


# Differential Expression Analysis Workflow

After running your Quantification software, DIA-NN, MAXQUANT, FragPipe-TMT, FragPipe-DIA or FragPipe-LFQ,
the quantification results are in an `data_dir`.
Please add the `.fasta` file which was used by the quantification software to the `data_dir`.

## Workflow Overview

1. [Create Dataset](#1-create-dataset)
2. [Generate Quality Control (QC)](#2-generate-quality-control-qc)
3. [Generate prolfqua YAML](#3-generate-prolfqua-yaml)
4. [Run Differential Expression Analysis](#4-run-differential-expression-analysis)


## 1. Create Dataset

The first step involves preparing the dataset by providing the experiment annotation. This is done using the `prolfqua_dataset.sh` script.

- **Input**: directory containing identification/quantification software outputs
- **Output**: csv, tsv or xlsx file template,

To create a `prolfquapp` compatible experiment annotation file run:

```bash
./prolfqua_dataset.sh -i data_dir/ -s DIANN -d annotation.xlsx
```

The `annotation.xlsx` file will be generated, and will contain 5 columns. 

- Relative.Path/Path/raw.file/channel/ (unique*)
- name - used in tables and figures (unique*)
- group/experiment/ - main factor
- subject/bioreplicate (optional**) - blocking factor
- control - used to specify the control condition (C) (optional)

* The rows must contain a unique value (no duplicates per column)
** If the experiment is not paired, or has no blocking factor (e.g. batch, cell line) delete the subject column.

The column raw.file is already filled out, based on the information available in the input directory.
You will need to fill out the missing columns, e.g. group, subject, control. The column names are not case sensitive.


The `annotation.xlsx` file will be generated, and you will need to fill out the empty cells.

## 2. Generate Quality Control (QC)

The `prolfqua_qc.sh` script will create a QC report. The report consists of two HTML documents and XLSX file. 

- **Input**: Dataset from step 1 and directory containing identification/quantification software outputs
- **Output**: QC report and visualizations

```bash
./prolfqua_qc.sh -i data_dir/ -p ProjectName -O ordername -w WorkunitName -d annotation.xlsx -s DIANN -o where_to_write_results
```

This will generate a subfolder which starts with "QC_" with all the analysis results.

## 3. Generate prolfquapp yaml

Using the `./prolfqua_yaml.R` command line tool you can set the parameters of the DEA and create a configuration file in YAML format.

- **Output**: Yaml file 


```
./prolfqua_yaml.sh -y config.yaml
```

To see which parameters can be set using `prolfqua_yaml.sh` use the `-h` switch. Other parameters you can set by editing the yaml file.


## 4. Run Differential Expression Analysis

Finally, the `prolfqua_dea.sh` script runs the differential expression analysis using the configuration file generated in the previous step.

- **Input**: directory containing identification/quantification software outputs, dataset from step 1, configuration file from step 3.
- **Output**: folder containing differential expression analysis results (html files, excel tables, rank files, SummarizedExperiment.rds).

After setting the parameters in the config.yaml file you can run the DEA analysis by:

```
./prolfqua_dea.sh -i data_dir/ -d annotation.xlsx -y config.yaml -w NameOfAnalysis -s DIANN
```

This will generate a subfolder which starts with "DEA_" and writes all the analysis results as well as the input data.


## ASMS poster:  Streamlining Protein Differential Expression Analysis in Core Facilities

![prolfquapp_ASMS_poster](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/prolfquapp_PosterPNG.png?raw=true)

### Introduction

The prolfquapp is a user-friendly application to streamline protein differential expression analysis (DEA) in core facilities. The application leverages the preprocessing methods and statistical models implemented in the R package [prolfqua](https://github.com/fgcz/prolfqua). It generates dynamic HTML reports that contain quality control plots and visualizations of the DEA results. The prolfquapp also exports results in multiple formats, including XLSX files, .rnk or txt files for gene set enrichment analysis, and Bioconductor SummarizedExperiment format for import into interactive visualization tools like OmicsViewer and iSEE. The prolfquapp offers a comprehensive and efficient solution for researchers seeking to analyze protein differential expression data in a core facility setting.

### Methods

The prolfquapp is a highly configurable application. The application interfaces with the laboratory data management systems through a YAML file and a '.tsv' file. To showcase this functionality, we developed an R Shiny application that collects user inputs, generates the configuration files, and runs prolfquapp. The YAML file enables the specification of key parameters, such as data normalization and modeling methods. It can also contain information about the input data and project integrated into the HTML report. This makes prolfquapp a user-friendly and flexible tool for protein DEA.

### Preliminary data

We developed an application based on the *R* package for users unfamiliar with statistics and *R* programming to ease the usage barriers of the *prolfqua* package. Furthermore, we integrated it into the data management platform B-Fabric. This integration enables users to select the input data and basic settings in a graphical user interface (GUI).
Using sample annotation information stored in B-Fabric, the user creates a table with the following columns: the sample name, a column with the sample groupings, an optional column with an additional explanatory variable, and finally, a column that defines which group differences to compute. This file format is similar to the one used by SAINTexpress, MSstats, or LFQAnalyst. In addition, other parameters, e.g., which of the modeling methods implemented in the prolfqua R-package can be specified.
The application then generates multiple result files, including an HTML report. The report integrates project-related information from B-Fabric and introduces the DEA. Furthermore, QC plots and summaries focus on the quality of protein and peptide identification and quantification. Next, we show the results of the DEA using dynamic tables and figures, which can be queried and filtered. Finally, in the section on additional analysis, we introduce methods of downstream omics analysis, e.g., gene set enrichment analysis, and provide links to web services (e.g., Webgestalt, String-db). 
Most importantly, the user receives all the data and code to reproduce the analysis on his infrastructure. This way, prolfquapp, and B-Fabric help scientists meet requirements from funding agencies, journals, and academic institutions while publishing their data according to the FAIR data principles. The source code of the prolfquapp R package is available from https://github.com/wolski/prolfquapp.

### Novel aspect

The prolfquapp can be tightly integrated with the LIMS system, and, at the same time, You can replicate the analysis on your laptop computer. 


# Implementation

- Yaml file is generated either using a Shiny application or a command line application.
- The parameters of the application are stored in an R list

## Yaml file



## Annotation file

A data frame with the sample annotation must be provided

- column which allows to map channel or raw file names to the sample name and explanatory variables.
  - if TMT analysis workflow one column starting with "channel" (either upper or lower case), with the channel (FragPipeTMT)
  - if DIA or DDA workflow then one column named exactly "Relative.Path", containing the raw filenames. Eeach  file name must be unique.
- one column starting with "^name" which is a unique sample name. These names will be used to label figures.
- one column starting with "group|^bait|^Experiment" (either upper or lower case). 
- Optional : no or exactly one column name starting with "^subject|^BioReplicate" (either upper or lower case)
- Optional : no or exactly one column name starting with "Contrast" (either upper or lower case)

## Raw input files


# Related software

- Einprot https://github.com/fmicompbio/einprot
- LFQAnalyst https://analyst-suite.monash-proteomics.cloud.edu.au/apps/lfq-analyst/ and https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst
- POMAShiny https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009148
- MSDap https://github.com/ftwkoopmans/msdap
