# Create DEA report in html and write data to xlsx table

For use examples see run_scripts directory

## Usage

``` r
make_DEA_report(lfqdata, protAnnot, GRP2)
```

## Arguments

- lfqdata:

  LFQData object

- protAnnot:

  ProteinAnnotation object

- GRP2:

  list with named arguments i.e. Contrasts, projectID, projectName,
  workunitID, nrPeptides, Diffthreshold, FDRthreshold

## See also

Other workflow:
[`make_SummarizedExperiment()`](https://prolfqua.github.io/prolfquapp/reference/make_SummarizedExperiment.md)
