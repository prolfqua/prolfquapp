# preprocess_dummy

preprocess_dummy

## Usage

``` r
preprocess_dummy(
  quant_data,
  fasta_file,
  annotation,
  pattern_contaminants = "^zz|^CON|Cont_",
  pattern_decoys = "^REV_|^rev",
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

- nr_peptides:

  accepted for interface uniformity but ignored (interface stub)
