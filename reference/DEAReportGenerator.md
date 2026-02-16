# DEAReportGenerator

DEAReportGenerator

DEAReportGenerator

## Public fields

- `lfqdata`:

  LFQData object containing the quantitative data

- `GRP2`:

  ProlfquAppConfig object containing analysis configuration

- `prot_annot`:

  ProteinAnnotation object

- `Contrasts`:

  list of contrasts for differential expression analysis

- `fname`:

  filename for DEA results

- `qcname`:

  filename for QC results

- `resultdir`:

  directory for storing results

- `ZIPDIR`:

  zip directory path

## Methods

### Public methods

- [`DEAReportGenerator$new()`](#method-DEAReportGenerator-new)

- [`DEAReportGenerator$write_DEA_all()`](#method-DEAReportGenerator-write_DEA_all)

- [`DEAReportGenerator$render_DEA()`](#method-DEAReportGenerator-render_DEA)

- [`DEAReportGenerator$make_boxplots()`](#method-DEAReportGenerator-make_boxplots)

- [`DEAReportGenerator$prep_result_list()`](#method-DEAReportGenerator-prep_result_list)

- [`DEAReportGenerator$make_SummarizedExperiment()`](#method-DEAReportGenerator-make_SummarizedExperiment)

- [`DEAReportGenerator$clone()`](#method-DEAReportGenerator-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize DEAReportGenerator with data and configuration

#### Usage

    DEAReportGenerator$new(lfqdata, GRP2, prot_annot, Contrasts)

#### Arguments

- `lfqdata`:

  LFQData object containing quantitative data

- `GRP2`:

  ProlfquAppConfig object with analysis configuration

- `prot_annot`:

  ProteinAnnotation object

- `Contrasts`:

  list of contrasts for differential expression analysis

------------------------------------------------------------------------

### Method [`write_DEA_all()`](https://prolfqua.github.io/prolfquapp/reference/write_DEA_all.md)

Write all DEA results to files

#### Usage

    DEAReportGenerator$write_DEA_all()

------------------------------------------------------------------------

### Method [`render_DEA()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_report.md)

Render DEA report using R Markdown

#### Usage

    DEAReportGenerator$render_DEA(
      htmlname,
      markdown = "_Grp2Analysis.Rmd",
      word = FALSE
    )

#### Arguments

- `htmlname`:

  name for the output HTML file

- `markdown`:

  path to the R Markdown template file

- `word`:

  logical, if TRUE output Word document, otherwise HTML

------------------------------------------------------------------------

### Method `make_boxplots()`

Generate boxplots for quality control

#### Usage

    DEAReportGenerator$make_boxplots()

------------------------------------------------------------------------

### Method `prep_result_list()`

Prepare result list with all analysis outputs

#### Usage

    DEAReportGenerator$prep_result_list()

#### Returns

list containing all analysis results

------------------------------------------------------------------------

### Method [`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md)

Create SummarizedExperiment object from analysis results

#### Usage

    DEAReportGenerator$make_SummarizedExperiment(
      strip = "~lfq~light",
      .url_builder = bfabric_url_builder
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
