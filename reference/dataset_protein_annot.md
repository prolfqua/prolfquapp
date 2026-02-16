# Dataset protein annot

Dataset protein annot

## Usage

``` r
dataset_protein_annot(
  msdata,
  idcol = c(protein_Id = "Protein.Group"),
  protein_annot = "fasta.header",
  more_columns = c("nrPeptides", "fasta.id")
)
```

## Arguments

- msdata:

  data frame

- idcol:

  named vector mapping protein ID column

- protein_annot:

  fasta header column name

- more_columns:

  more columns to include
