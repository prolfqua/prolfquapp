# Make duplicated CompoundDiscoverer feature/sample rows unique

Make duplicated CompoundDiscoverer feature/sample rows unique

## Usage

``` r
make_CD_duplicate_features_unique(
  data,
  feature_col = "Feature_ID",
  sample_col = "Sample"
)
```

## Arguments

- data:

  long CompoundDiscoverer data frame

- feature_col:

  feature identifier column

- sample_col:

  sample identifier column

## Value

data frame with unique feature identifiers for affected rows
