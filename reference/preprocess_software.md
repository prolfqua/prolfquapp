# collects preprocess methods for various software

collects preprocess methods for various software

## Usage

``` r
preprocess_software(
  indir,
  annotation,
  preprocess_functions,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^rev_|^REV_",
  nr_peptides = 1,
  extreader = NULL
)
```

## Arguments

- indir:

  input directory with quantification data

- annotation:

  annotation list from read_annotation

- preprocess_functions:

  list or R6 object with get_files, preprocess, extra_args

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

- nr_peptides:

  minimum number of distinct (stripped) peptides per parent protein;
  forwarded only to readers that declare an \`nr_peptides\` argument.
  Readers without it are warned about (and left unfiltered) when
  \`nr_peptides \> 1\`. Default 1 (no filtering).

- extreader:

  optional external reader configuration

## Examples

``` r
# example code
annot <- data.frame(
  file = c("a1.raw", "a2.raw", "a3.raw", "a4.raw"),
  name = c("aa", "ba", "aa", "ba"),
  group = c("a", "a", "b", "b")
)

annot <- read_annotation(annot, QC = TRUE)
#> INFO [2026-07-02 13:50:29] Using derived sample display names in column 'sampleName'.
preprocess_functions <- prolfquapp::prolfqua_preprocess_functions[["DUMMY"]]
res <- preprocess_software(".", annot, preprocess_functions)
#> INFO [2026-07-02 13:50:29] Files data: data.path
#> INFO [2026-07-02 13:50:29] Files fasta: fasta.files.path

xx <- prolfquapp::ExternalReader$new()
xx$extra_args <- "list()"
xx$get_files <- "prolfquapp::get_dummy_files"
xx$preprocess <- "prolfquapp::preprocess_dummy"
res <- preprocess_software(".", annotation = annot, preprocess_functions = xx)
#> INFO [2026-07-02 13:50:29] Files data: data.path
#> INFO [2026-07-02 13:50:29] Files fasta: fasta.files.path
xx <- prolfquapp::ExternalReader$new()
```
