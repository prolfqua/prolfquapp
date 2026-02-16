# preprocess MQ peptide

preprocess MQ peptide

preprocess MQ peptide

## Usage

``` r
preprocess_MQ_peptide(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1
)

preprocess_MQ_peptide(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev_",
  hierarchy_depth = 1
)
```

## Arguments

- quant_data:

  path to peptides.txt file

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
