# preprocess DIANN ouput, filter by q_value and nr_peptides

preprocess DIANN ouput, filter by q_value and nr_peptides

preprocess DIANN ouput, filter by q_value and nr_peptides

## Usage

``` r
preprocess_BGS(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 2
)

preprocess_BGS(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
  q_value = 0.01,
  hierarchy_depth = 2
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

## Value

list with lfqdata and protein annotation

list with lfqdata and protein annotation

## Examples

``` r
if (FALSE) { # \dontrun{
x <- get_BGS_files("DefaultParsing")
bgs <- read_BGS(x$data)
annot <- data.frame(raw.file = bgs$R.FileName |> unique(),
 Name = paste(c(rep("A",3),rep("B",3)),1:6, sep="_"),
group = c(rep("A",3),rep("B",3)))
annotation <- annot |> prolfquapp::read_annotation(QC = TRUE)
#debug(preprocess_BGS)
xd <- preprocess_BGS(x$data, x$fasta, annotation)
} # }
if (FALSE) { # \dontrun{
x <- get_BGS_files("DefaultParsing")
bgs <- read_BGS(x$data)
annot <- data.frame(raw.file = bgs$R.FileName |> unique(),
 Name = paste(c(rep("A",3),rep("B",3)),1:6, sep="_"),
group = c(rep("A",3),rep("B",3)))
annotation <- annot |> prolfquapp::read_annotation(QC = TRUE)
#debug(preprocess_BGS)
xd <- preprocess_BGS(x$data, x$fasta, annotation)
} # }
```
