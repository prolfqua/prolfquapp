# preprocess DIANN ouput, filter by q_value and nr_peptides

preprocess DIANN ouput, filter by q_value and nr_peptides

## Usage

``` r
preprocess_DIANN(
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

list with lfqdata and protein annotation

## Examples

``` r
if (FALSE) { # \dontrun{
x <- get_DIANN_files("inst/application/DIANN/2706527/")

annotation <- file.path("inst/application/DIANN/2706527/dataset.csv") |>
  readr::read_csv() |>
  prolfquapp::read_annotation(QC = TRUE)
x$fasta
undebug(preprocess_DIANN)
xd <- preprocess_DIANN(x$data, x$fasta, annotation)
xd$lfqdata$hierarchy_counts()
xd <- preprocess_DIANN(x$data, x$fasta, annotation, nr_peptides = 2)
xd$lfqdata$hierarchy_counts()
} # }
```
