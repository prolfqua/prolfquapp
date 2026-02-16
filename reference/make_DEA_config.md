# create GRP2 configuration. Use this function if there is no Yaml Input.

create GRP2 configuration. Use this function if there is no Yaml Input.

## Usage

``` r
make_DEA_config(
  ZIPDIR = ".",
  PROJECTID = "",
  ORDERID = "",
  WORKUNITID = "",
  Normalization = c("none", "vsn", "quantile", "robscale"),
  aggregation = c("medpolish", "top3", "lmrob"),
  Diffthreshold = 1,
  FDRthreshold = 0.1,
  removeContaminants = FALSE,
  removeDecoys = FALSE,
  patternDecoys = "^REV_",
  patternContaminants = "^zz",
  application = "FragPipeTMT",
  nrPeptides = 2
)
```

## Arguments

- ZIPDIR:

  output zip directory

- PROJECTID:

  project identifier

- ORDERID:

  order identifier

- WORKUNITID:

  workunit identifier

- Normalization:

  normalization method

- aggregation:

  aggregation method

- Diffthreshold:

  fold-change difference threshold

- FDRthreshold:

  FDR significance threshold

- removeContaminants:

  if TRUE remove contaminants

- removeDecoys:

  if TRUE remove decoy sequences

- patternDecoys:

  default "^REV\_"

- patternContaminants:

  default "^zz\_"

- application:

  software application name

- nrPeptides:

  minimum number of peptides per protein

## Examples

``` r
DEAconfig <- make_DEA_config()
#> Warning: DEPRECATED
DEAconfig
#> $project_spec
#> $project_spec$project_Id
#> [1] ""
#> 
#> $project_spec$project_name
#> [1] ""
#> 
#> $project_spec$order_Id
#> [1] ""
#> 
#> $project_spec$workunit_Id
#> [1] ""
#> 
#> 
#> $pop
#> $pop$transform
#> [1] "none"
#> 
#> $pop$aggregate
#> [1] "medpolish"
#> 
#> $pop$Diffthreshold
#> [1] 1
#> 
#> $pop$FDRthreshold
#> [1] 0.1
#> 
#> $pop$removeCon
#> [1] FALSE
#> 
#> $pop$removeDecoys
#> [1] FALSE
#> 
#> $pop$revpattern
#> [1] "^REV_"
#> 
#> $pop$contpattern
#> [1] "^zz"
#> 
#> $pop$nr_peptdes
#> [1] 2
#> 
#> 
#> $Software
#> [1] "FragPipeTMT"
#> 
#> $zipdir
#> [1] "."
#> 
```
