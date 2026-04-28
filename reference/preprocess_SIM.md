# Preprocess simulated data

Returns a simulated LFQData + ProteinAnnotation, bypassing all file I/O.
Designed for integration testing of CMD scripts via `--software SIM`.

## Usage

``` r
preprocess_SIM(
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

  ignored (placeholder)

- fasta_file:

  ignored (placeholder)

- annotation:

  annotation list from
  [`read_annotation`](https://prolfqua.github.io/prolfquapp/reference/read_annotation.md)

- pattern_contaminants:

  regex for contaminant proteins

- pattern_decoys:

  regex for decoy proteins

- hierarchy_depth:

  1 = protein level, 2 = peptide level

## Value

list with `lfqdata` (LFQData) and `protein_annotation`
(ProteinAnnotation)

## Details

The simulated data is reconfigured to use the annotation's factor prefix
so that contrasts derived from the annotation match the model terms.
