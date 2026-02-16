# R6 class representing ProlfquApp configuration

R6 class representing ProlfquApp configuration

R6 class representing ProlfquApp configuration

## See also

Other ProlfquAppConfig:
[`ExternalReader`](https://prolfqua.github.io/prolfquapp/reference/ExternalReader.md),
[`ProcessingOptions`](https://prolfqua.github.io/prolfquapp/reference/ProcessingOptions.md),
[`ProjectSpec`](https://prolfqua.github.io/prolfquapp/reference/ProjectSpec.md),
[`make_DEA_config_R6()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_config_R6.md),
[`set_list_to_R6()`](https://prolfqua.github.io/prolfquapp/reference/set_list_to_R6.md)

## Public fields

- `processing_options`:

  ProcessingOption R6 class

- `project_spec`:

  Project Spec R6 class

- `software`:

  name of input software

- `prefix`:

  either QC or DEA

- `zipdir_name`:

  results should go to zipdir_name

- `path`:

  path to working directory

- `pop`:

  optional processing options

- `RES`:

  results

- `group`:

  group prefix

- `ext_reader`:

  external reader configuration

## Methods

### Public methods

- [`ProlfquAppConfig$new()`](#method-ProlfquAppConfig-new)

- [`ProlfquAppConfig$set_zipdir_name()`](#method-ProlfquAppConfig-set_zipdir_name)

- [`ProlfquAppConfig$get_zipdir()`](#method-ProlfquAppConfig-get_zipdir)

- [`ProlfquAppConfig$get_result_dir()`](#method-ProlfquAppConfig-get_result_dir)

- [`ProlfquAppConfig$get_input_dir()`](#method-ProlfquAppConfig-get_input_dir)

- [`ProlfquAppConfig$as_list()`](#method-ProlfquAppConfig-as_list)

- [`ProlfquAppConfig$clone()`](#method-ProlfquAppConfig-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize ProlfquAppConfig with processing options, project spec, and
external reader

#### Usage

    ProlfquAppConfig$new(
      processing_options,
      project_spec,
      ext_reader,
      zipdir_name = ".",
      path = ".",
      software = "DIANN",
      prefix = "DEA"
    )

#### Arguments

- `processing_options`:

  instance of ProcessingOptions

- `project_spec`:

  instance of ProjectSpec

- `ext_reader`:

  instance of ExternalReader

- `zipdir_name`:

  where to store results

- `path`:

  working directory path

- `software`:

  name of input software

- `prefix`:

  either QC or DEA

------------------------------------------------------------------------

### Method `set_zipdir_name()`

Set zip directory name based on project information and date

#### Usage

    ProlfquAppConfig$set_zipdir_name()

#### Returns

the generated zip directory name

------------------------------------------------------------------------

### Method `get_zipdir()`

Get the full path to the zip directory

#### Usage

    ProlfquAppConfig$get_zipdir()

#### Returns

full path to zip directory

------------------------------------------------------------------------

### Method `get_result_dir()`

Get the results directory path

#### Usage

    ProlfquAppConfig$get_result_dir()

#### Returns

path to results directory

------------------------------------------------------------------------

### Method `get_input_dir()`

Get the input directory path

#### Usage

    ProlfquAppConfig$get_input_dir()

#### Returns

path to input directory

------------------------------------------------------------------------

### Method `as_list()`

Convert R6 object to list

#### Usage

    ProlfquAppConfig$as_list()

#### Returns

list representation of the R6 object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ProlfquAppConfig$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

r6obj_config <- ProlfquAppConfig$new(ProcessingOptions$new(), ProjectSpec$new(), ExternalReader$new())
xx <- prolfqua::R6_extract_values(r6obj_config)
yaml::write_yaml(xx, file = file.path(tempdir(), "test.yaml"))
config <- yaml::read_yaml(file = file.path(tempdir(), "test.yaml"))

r6obj_config$set_zipdir_name()
#> [1] "DEA_20260216_vsn"

r6obj_config$get_zipdir()
#> [1] "./DEA_20260216_vsn"
r6obj_config$get_result_dir()
#> [1] "./DEA_20260216_vsn/Results_WU_"
r6obj_config$get_input_dir()
#> [1] "./DEA_20260216_vsn/Inputs_WU_"
```
