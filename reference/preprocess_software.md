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
preprocess_functions <- prolfquapp::prolfqua_preprocess_functions[["DUMMY"]]
res <- preprocess_software(".", annot, preprocess_functions)
#> INFO [2026-02-25 16:36:05] Files data: data.path
#> INFO [2026-02-25 16:36:05] Files fasta: fasta.files.path

xx <- prolfquapp::ExternalReader$new()
xx$extra_args <- "list()"
xx$get_files <- "prolfquapp::get_dummy_files"
xx$preprocess <- "prolfquapp::preprocess_dummy"
res <- preprocess_software(".", annotation = annot, preprocess_functions = xx)
#> INFO [2026-02-25 16:36:05] Files data: data.path
#> INFO [2026-02-25 16:36:05] Files fasta: fasta.files.path
xx <- prolfquapp::ExternalReader$new()
```
