# preprocess mzMine input

preprocess mzMine input

## Usage

``` r
preprocess_mzMine(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = NULL,
  pattern_decoys = NULL,
  annotated = FALSE
)
```

## Arguments

- quant_data:

  path to mzMine features csv file

- fasta_file:

  path to annotations csv file

- annotation:

  annotation list from read_annotation

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

- annotated:

  if TRUE only keep annotated features

## Examples

``` r
if(FALSE){
xd <- "outputs-20250407T1707/bfabric/input_dataset.tsv"
annot <- readr::read_tsv(xd)

annotation <- read_annotation(annot, QC = TRUE)
xd <- "outputs-20250407T1707/"
files <- get_mzMine_files(path)
files
undebug(preprocess_mzMine)
res <- preprocess_mzMine(files$data, files$fasta , annotation)
dim(res$lfqdata$data)
res <- preprocess_mzMine(files$data, files$fasta , annotation, annotated = TRUE)
dim(res$lfqdata$data)
}
```
