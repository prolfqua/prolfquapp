# build Dataset protein annot, defaults are compatible with DIANN

build Dataset protein annot, defaults are compatible with DIANN

## Usage

``` r
build_protein_annot(
  lfqdata,
  msdata,
  idcol = c(protein_Id = "Protein.Group"),
  cleaned_protein_id = "Protein.Group.2",
  protein_description = "fasta.header",
  exp_nr_children = "nrPeptides",
  full_id = "fasta.id",
  more_columns = c("fasta.id"),
  pattern_contaminants = "^zz|^CON",
  pattern_decoys = "REV_"
)
```

## Arguments

- lfqdata:

  LFQData

- msdata:

  data frame

- idcol:

  named vector mapping protein ID column

- cleaned_protein_id:

  column with cleaned protein ID

- protein_description:

  column with protein description

- exp_nr_children:

  column with number of peptides

- full_id:

  column with full protein ID

- more_columns:

  additional columns to include

- pattern_contaminants:

  regex pattern for contaminants

- pattern_decoys:

  regex pattern for decoys

## Examples

``` r
# example code
```
