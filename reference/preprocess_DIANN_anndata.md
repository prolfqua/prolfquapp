# Preprocess DIANN output and return AnnData

Same interface as
[`preprocess_DIANN`](https://prolfqua.github.io/prolfquapp/reference/preprocess_DIANN.md)
but returns an
[`anndataR::AnnData`](https://anndataR.scverse.org/reference/AnnData.html)
object instead of `list(lfqdata, protein_annotation)`.

## Usage

``` r
preprocess_DIANN_anndata(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 1,
  nr_peptides = 1
)
```

## Arguments

- quant_data:

  path to quantification data file

- fasta_file:

  path to fasta file(s)

- annotation:

  annotation list from read_annotation

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

- q_value:

  q-value threshold for filtering

- hierarchy_depth:

  hierarchy depth for aggregation

- nr_peptides:

  minimum number of peptides per protein

## Value

[`anndataR::AnnData`](https://anndataR.scverse.org/reference/AnnData.html)
object

## Details

The AnnData `uns` slot contains three namespaces:

- X_layer_name:

  Name of the primary intensity column stored in X

- exploreDE:

  Column role metadata compatible with anndata_omics_bridge

- prolfquapp:

  Round-trip reconstruction metadata (config, protein annotation)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- get_DIANN_files("inst/application/DIANN/2706527/")
annotation <- file.path("inst/application/DIANN/2706527/dataset.csv") |>
  readr::read_csv() |>
  prolfquapp::read_annotation(QC = TRUE)
adata <- preprocess_DIANN_anndata(x$data, x$fasta, annotation)
} # }
```
