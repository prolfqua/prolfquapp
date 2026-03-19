# Convert LFQData + ProteinAnnotation to AnnData with round-trip metadata

Builds an AnnData object from LFQData and ProteinAnnotation, with the
`prolfquapp` uns namespace for lossless round-trip via
[`LFQData_from_anndata`](https://prolfqua.github.io/prolfquapp/reference/LFQData_from_anndata.md).

## Usage

``` r
preprocess_anndata_from_lfq(lfqdata, protAnnot, source_software = "unknown")
```

## Arguments

- lfqdata:

  LFQData object

- protAnnot:

  ProteinAnnotation object

- source_software:

  character, name of the source software (e.g. "DIANN")

## Value

[`anndataR::AnnData`](https://anndataR.scverse.org/reference/AnnData.html)
object

## Details

Handles both protein-level (hierarchyDepth=1) and peptide-level
(hierarchyDepth=2) data. For peptide-level data, var contains one row
per feature (peptide) with protein annotation joined in.

## Examples

``` r
if (requireNamespace("anndataR", quietly = TRUE)) {
  # Protein-level
  res <- sim_data_protAnnot(Nprot = 10, PROTEIN = TRUE)
  adata <- preprocess_anndata_from_lfq(res$lfqdata, res$pannot, "simulated")
  adata  # 12 obs x 10 var

  # Peptide-level (multiple peptides per protein)
  res2 <- sim_data_protAnnot(Nprot = 10, PROTEIN = FALSE)
  adata2 <- preprocess_anndata_from_lfq(res2$lfqdata, res2$pannot, "simulated")
  adata2  # 12 obs x ~28 var (peptides)

  # Round-trip back to LFQData
  back <- LFQData_from_anndata(adata2)
  back$lfqdata$hierarchy_counts()
}
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function
#> # A tibble: 1 × 3
#>   isotopeLabel protein_Id peptide_Id
#>   <chr>             <int>      <int>
#> 1 light                10         28
```
