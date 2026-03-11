# DEAReportGenerator

DEAReportGenerator

DEAReportGenerator

## Details

Generates all output files for a differential expression analysis. Uses
a DEAnalyse object as data source instead of the legacy GRP2\$RES list.

## Public fields

- `deanalyse`:

  DEAnalyse object containing all analysis results

- `GRP2`:

  ProlfquAppConfig object containing analysis configuration

- `fname`:

  filename prefix for DEA results

- `qcname`:

  filename prefix for QC results

- `resultdir`:

  directory for storing results

- `ZIPDIR`:

  zip directory path

## Methods

### Public methods

- [`DEAReportGenerator$new()`](#method-DEAReportGenerator-new)

- [`DEAReportGenerator$prep_result_list()`](#method-DEAReportGenerator-prep_result_list)

- [`DEAReportGenerator$write_DEA()`](#method-DEAReportGenerator-write_DEA)

- [`DEAReportGenerator$render_DEA()`](#method-DEAReportGenerator-render_DEA)

- [`DEAReportGenerator$make_boxplots()`](#method-DEAReportGenerator-make_boxplots)

- [`DEAReportGenerator$filter_data()`](#method-DEAReportGenerator-filter_data)

- [`DEAReportGenerator$get_protein_boxplots()`](#method-DEAReportGenerator-get_protein_boxplots)

- [`DEAReportGenerator$contrasts_to_Grob()`](#method-DEAReportGenerator-contrasts_to_Grob)

- [`DEAReportGenerator$get_protein_boxplots_contrasts()`](#method-DEAReportGenerator-get_protein_boxplots_contrasts)

- [`DEAReportGenerator$write_protein_boxplots()`](#method-DEAReportGenerator-write_protein_boxplots)

- [`DEAReportGenerator$write_DEA_all()`](#method-DEAReportGenerator-write_DEA_all)

- [`DEAReportGenerator$make_SummarizedExperiment()`](#method-DEAReportGenerator-make_SummarizedExperiment)

- [`DEAReportGenerator$clone()`](#method-DEAReportGenerator-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize DEAReportGenerator

#### Usage

    DEAReportGenerator$new(deanalyse, GRP2, name = "")

#### Arguments

- `deanalyse`:

  DEAnalyse R6 object with completed analysis

- `GRP2`:

  ProlfquAppConfig R6 object

- `name`:

  optional name prefix for output files

------------------------------------------------------------------------

### Method `prep_result_list()`

Prepare result list with all analysis outputs for XLSX

#### Usage

    DEAReportGenerator$prep_result_list()

#### Returns

list containing all analysis results (14 sheets)

------------------------------------------------------------------------

### Method `write_DEA()`

Write DEA results (XLSX, ORA, GSEA files)

#### Usage

    DEAReportGenerator$write_DEA(ORA = TRUE, GSEA = TRUE)

#### Arguments

- `ORA`:

  if TRUE write ORA gene lists

- `GSEA`:

  if TRUE write GSEA rank files

#### Returns

list with xlsx_file, ora_files, gsea_files paths

------------------------------------------------------------------------

### Method `render_DEA()`

Render DEA report using R Markdown

#### Usage

    DEAReportGenerator$render_DEA(
      htmlname,
      markdown = "_Grp2Analysis_V2_R6.Rmd",
      word = FALSE,
      toc = TRUE
    )

#### Arguments

- `htmlname`:

  name for the output HTML file

- `markdown`:

  path to the R Markdown template file

- `word`:

  logical, if TRUE output Word document, otherwise HTML

- `toc`:

  logical, if TRUE include table of contents

#### Returns

path to the output file

------------------------------------------------------------------------

### Method `make_boxplots()`

Generate sample-level boxplots for quality control

#### Usage

    DEAReportGenerator$make_boxplots(boxplot = TRUE)

#### Arguments

- `boxplot`:

  logical, if TRUE write boxplots

------------------------------------------------------------------------

### Method `filter_data()`

Get subset of transformed data for significant proteins

#### Usage

    DEAReportGenerator$filter_data()

------------------------------------------------------------------------

### Method `get_protein_boxplots()`

Get per-protein boxplots for significant proteins

#### Usage

    DEAReportGenerator$get_protein_boxplots()

------------------------------------------------------------------------

### Method `contrasts_to_Grob()`

Convert significant contrast results to table grobs

#### Usage

    DEAReportGenerator$contrasts_to_Grob()

------------------------------------------------------------------------

### Method `get_protein_boxplots_contrasts()`

Get per-protein boxplots combined with contrast summary tables

#### Usage

    DEAReportGenerator$get_protein_boxplots_contrasts()

------------------------------------------------------------------------

### Method `write_protein_boxplots()`

Write per-protein boxplots with contrast tables to PDF

#### Usage

    DEAReportGenerator$write_protein_boxplots(filename = "boxplots")

#### Arguments

- `filename`:

  base filename (without extension)

------------------------------------------------------------------------

### Method `write_DEA_all()`

Write all DEA results: XLSX, ORA, GSEA, HTML reports, boxplots

#### Usage

    DEAReportGenerator$write_DEA_all(
      boxplot = TRUE,
      render = TRUE,
      ORA = TRUE,
      GSEA = TRUE,
      markdown = "_Grp2Analysis_V2_R6.Rmd",
      markdown_qc = "_DiffExpQC_R6.Rmd",
      toc = TRUE
    )

#### Arguments

- `boxplot`:

  if TRUE generate boxplots

- `render`:

  if TRUE render HTML reports

- `ORA`:

  if TRUE write ORA gene lists

- `GSEA`:

  if TRUE write GSEA rank files

- `markdown`:

  Rmd template for main DEA report

- `markdown_qc`:

  Rmd template for QC report

- `toc`:

  if TRUE include table of contents

#### Returns

list with dea_file, qc_file, data_files paths

------------------------------------------------------------------------

### Method [`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md)

Create SummarizedExperiment object from analysis results

#### Usage

    DEAReportGenerator$make_SummarizedExperiment(
      strip = "~lfq~light",
      .url_builder = prolfquapp::bfabric_url_builder
    )

#### Arguments

- `strip`:

  pattern to strip from rownames

- `.url_builder`:

  function to build URLs for bfabric

#### Returns

SummarizedExperiment object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DEAReportGenerator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
