# Convert AnnData back to LFQData + ProteinAnnotation

Reconstructs the `list(lfqdata, protein_annotation)` pair from an
AnnData object that was created by
[`preprocess_DIANN_anndata`](https://prolfqua.github.io/prolfquapp/reference/preprocess_DIANN_anndata.md)
(or any producer that writes the `prolfquapp` uns namespace).

## Usage

``` r
LFQData_from_anndata(adata)
```

## Arguments

- adata:

  an
  [`anndataR::AnnData`](https://anndataR.scverse.org/reference/AnnData.html)
  object with `uns$prolfquapp`

## Value

list with `lfqdata` (LFQData) and `protein_annotation`
(ProteinAnnotation)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- sim_data_protAnnot()
adata <- preprocess_DIANN_anndata_from_lfq(res$lfqdata, res$pannot)
back <- LFQData_from_anndata(adata)
} # }
```
