# preprocess MSstats fragpipe

preprocess MSstats fragpipe

## Usage

``` r
preprocess_MSstats_FPDIA(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "",
  pattern_decoys = "",
  hierarchy_depth = 1
)
```

## Arguments

- quant_data:

  path to MSstats csv file

- fasta_file:

  path to fasta file(s)

- annotation:

  annotation list from read_annotation

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

- hierarchy_depth:

  hierarchy depth for aggregation
