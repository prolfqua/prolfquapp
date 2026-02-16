# external reader R6 class for handling external data sources

external reader R6 class for handling external data sources

external reader R6 class for handling external data sources

## See also

Other ProlfquAppConfig:
[`ProcessingOptions`](https://prolfqua.github.io/prolfquapp/reference/ProcessingOptions.md),
[`ProjectSpec`](https://prolfqua.github.io/prolfquapp/reference/ProjectSpec.md),
[`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md),
[`make_DEA_config_R6()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_config_R6.md),
[`set_list_to_R6()`](https://prolfqua.github.io/prolfquapp/reference/set_list_to_R6.md)

## Public fields

- `get_files`:

  function name for getting files

- `preprocess`:

  function name for preprocessing

- `extra_args`:

  extra arguments as string

- `dataset`:

  function name for getting dataset

## Methods

### Public methods

- [`ExternalReader$clone()`](#method-ExternalReader-clone)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExternalReader$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ExternalReader$new()
#> <ExternalReader>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     dataset: 
#>     extra_args: list()
#>     get_files: 
#>     preprocess: 
```
