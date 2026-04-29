# Preprocess a CompoundDiscoverer prolfqua ZIP export

Preprocess a CompoundDiscoverer prolfqua ZIP export

## Usage

``` r
preprocess_CD_export(long_file, sample_file, config, subset_column = NULL)
```

## Arguments

- long_file:

  path to \*\_long.csv

- sample_file:

  path to \*\_prolfqua_samples.csv

- config:

  ProlfquAppConfig object

- subset_column:

  optional column in long_file marking features to keep

## Value

list with lfqdata, protein_annotation, and annotation
