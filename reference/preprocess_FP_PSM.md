# preprocess FP psm, filter by purity_threshold and PeptideProphetProb

preprocess FP psm, filter by purity_threshold and PeptideProphetProb

## Usage

``` r
preprocess_FP_PSM(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  purity_threshold = 0.5,
  PeptideProphetProb = 0.9,
  hierarchy_depth = 1,
  parse_fun = tidy_FragPipe_psm
)
```

## Arguments

- quant_data:

  path to quantification data file(s)

- fasta_file:

  path to fasta file(s)

- annotation:

  annotation list from read_annotation

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

- purity_threshold:

  purity threshold for filtering

- PeptideProphetProb:

  PeptideProphet probability threshold

- hierarchy_depth:

  hierarchy depth for aggregation

- parse_fun:

  function for parsing PSM files

## Value

list with lfqdata and protein annotation
