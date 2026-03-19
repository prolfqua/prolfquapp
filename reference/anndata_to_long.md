# Convert AnnData wide matrix to long-format tibble

Internal helper. Melts the X matrix and any additional layers, joins
with obs (sample factors) and var (hierarchy columns).

## Usage

``` r
anndata_to_long(adata, config)
```

## Arguments

- adata:

  AnnData object

- config:

  AnalysisConfiguration

## Value

tibble in long format
