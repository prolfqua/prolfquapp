# create GRP2 configuration for differential expression analysis Use this function if there is no Yaml Input.

create GRP2 configuration for differential expression analysis Use this
function if there is no Yaml Input.

## Usage

``` r
make_DEA_config_R6(
  PATH = ".",
  PROJECTID = "",
  ORDERID = "",
  WORKUNITID = "",
  Normalization = c("none", "vsn", "quantile", "robscale"),
  aggregation = c("medpolish", "top3", "lmrob"),
  diff_threshold = 1,
  FDR_threshold = 0.1,
  nr_peptides = 1,
  removeContaminants = FALSE,
  removeDecoys = FALSE,
  patternDecoys = "^REV_|^rev_",
  patternContaminants = "^zz|^CON|Cont_",
  application = "DIANN",
  prefix = "DEA"
)
```

## Arguments

- PATH:

  working directory path

- PROJECTID:

  project identifier

- ORDERID:

  order identifier

- WORKUNITID:

  workunit identifier

- Normalization:

  normalization method: "none", "vsn", "quantile", "robscale"

- aggregation:

  aggregation method: "medpolish", "top3", "lmrob"

- diff_threshold:

  difference threshold

- FDR_threshold:

  FDR threshold

- nr_peptides:

  number of peptides required

- removeContaminants:

  should contaminants be removed

- removeDecoys:

  should decoys be removed

- patternDecoys:

  pattern for decoy proteins

- patternContaminants:

  pattern for contaminant proteins

- application:

  software application name

- prefix:

  analysis prefix (DEA or QC)

## Value

ProlfquAppConfig R6 object

## See also

Other ProlfquAppConfig:
[`ExternalReader`](https://prolfqua.github.io/prolfquapp/reference/ExternalReader.md),
[`ProcessingOptions`](https://prolfqua.github.io/prolfquapp/reference/ProcessingOptions.md),
[`ProjectSpec`](https://prolfqua.github.io/prolfquapp/reference/ProjectSpec.md),
[`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md),
[`set_list_to_R6()`](https://prolfqua.github.io/prolfquapp/reference/set_list_to_R6.md)

## Examples

``` r
DEAconfig <- make_DEA_config_R6(ORDERID = "1234", WORKUNITID = "1234")
DEAconfig$set_zipdir_name()
#> [1] "DEA_20260223_O1234_WU1234_none"
DEAconfig$get_zipdir()
#> [1] "./DEA_20260223_O1234_WU1234_none"
DEAconfig$get_result_dir()
#> [1] "./DEA_20260223_O1234_WU1234_none/Results_WU_1234"
DEAconfig$get_input_dir()
#> [1] "./DEA_20260223_O1234_WU1234_none/Inputs_WU_1234"
R6list <- prolfqua::R6_extract_values(DEAconfig)
```
