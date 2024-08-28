# prolfquapp: Generating Dynamic DEA Reports using a command line interface to the prolfqua R Package

Welcome to *prolfquapp* on GitHub! Here, you'll find everything you need to elevate your protein differential expression analysis.
Prolfquapp integrates powerful preprocessing methods and advanced statistical models from the prolfqua R package [prolfqua](https://github.com/fgcz/prolfqua)
[doi.org/10.1021/acs.jproteome.2c00441](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00441) to deliver insightful, clear visualizations and robust data outputs. 
Generate dynamic HTML reports, versatile file formats, and dive into interactive data visualization with ExploreDE.
Prolfquapp implements a command line interface to run protein differential expression analysis, that can be integrated into your workflow manager.

![prolfquapp](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/Prolfqapp_Highlight.png?raw=true)


# Differential Expression Analysis Workflow with prolfquapp

After running your Quantification software, DIA-NN, MAXQUANT, FragPipe-TMT, FragPipe-DIA or FragPipe-LFQ,
the quantification results are in an `data_dir`.
Please add the `.fasta` file which was used by the quantification software to the `data_dir`.

prolfquapp is a set of command line tool. To use it open you shell (linux, mac), or command window (windows).
Change into the directory with the identification/quantification results coming from FragPipe, MaxQuant, DIA-NN, Spectronaut etc.
In the directory execute:

```bash
R --vanilla -e "prolfquapp::copy_shell_script(workdir = '.')"
```

This will place the following four shell script files (linux), or bat files (windows) into your working directory:


```bash
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
- subject/bioreplicate (optional** or keep cells empty) - blocking factor
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

## How to install

**Linux**

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

**Docker**

## ASMS poster:  Streamlining Protein Differential Expression Analysis in Core Facilities

![prolfquapp_ASMS_poster](https://github.com/prolfqua/prolfquapp/blob/master/inst/poster/prolfquapp_PosterPNG.png?raw=true)


# Related software

- Einprot https://github.com/fmicompbio/einprot
- LFQAnalyst https://analyst-suite.monash-proteomics.cloud.edu.au/apps/lfq-analyst/ and https://github.com/MonashBioinformaticsPlatform/LFQ-Analyst
- POMAShiny https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009148
- MSDap https://github.com/ftwkoopmans/msdap
