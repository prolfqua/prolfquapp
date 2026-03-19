# Dataset protein annot

Extracts protein annotation from a data frame, renaming columns and
auto-detecting UniProt identifiers. For new code prefer
[`build_protein_annot`](https://prolfqua.github.io/prolfquapp/reference/build_protein_annot.md)
which returns a
[`ProteinAnnotation`](https://prolfqua.github.io/prolfquapp/reference/ProteinAnnotation.md)
R6 object.

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
