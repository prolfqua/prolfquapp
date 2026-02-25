# keep track of folder paths and create them if needed

keep track of folder paths and create them if needed

keep track of folder paths and create them if needed

## Public fields

- `outpath`:

  path

- `qc_dir`:

  qc results directory name

- `modelling_dir`:

  modeling results directory name

- `project_Id`:

  project_Id

- `order_Id`:

  order_Id

- `workunit_Id`:

  workunit_Id

- `inputData`:

  inputFile

- `inputAnnotation`:

  inputAnnotation xlsx

## Methods

### Public methods

- [`ProjectStructure$new()`](#method-ProjectStructure-new)

- [`ProjectStructure$create_outpath()`](#method-ProjectStructure-create_outpath)

- [`ProjectStructure$qc_path()`](#method-ProjectStructure-qc_path)

- [`ProjectStructure$modelling_path()`](#method-ProjectStructure-modelling_path)

- [`ProjectStructure$create()`](#method-ProjectStructure-create)

- [`ProjectStructure$reset()`](#method-ProjectStructure-reset)

- [`ProjectStructure$clone()`](#method-ProjectStructure-clone)

------------------------------------------------------------------------

### Method `new()`

create ProjectStructure

#### Usage

    ProjectStructure$new(
      outpath,
      project_Id,
      order_Id,
      workunit_Id,
      inputAnnotation,
      inputData,
      qc_dir = "qc_results",
      modelling_dir = "modelling_results"
    )

#### Arguments

- `outpath`:

  directory

- `project_Id`:

  bfabric project ID

- `order_Id`:

  bfabric order_Id

- `workunit_Id`:

  bfabric workunit_Id

- `inputAnnotation`:

  input annotation path

- `inputData`:

  input data path

- `qc_dir`:

  qc folder

- `modelling_dir`:

  modelling results folder

------------------------------------------------------------------------

### Method `create_outpath()`

create outpath

#### Usage

    ProjectStructure$create_outpath()

------------------------------------------------------------------------

### Method `qc_path()`

create qc dir

#### Usage

    ProjectStructure$qc_path(qc_dir)

#### Arguments

- `qc_dir`:

  QC directory

------------------------------------------------------------------------

### Method `modelling_path()`

create modelling path

#### Usage

    ProjectStructure$modelling_path(modelling_dir)

#### Arguments

- `modelling_dir`:

  directory with modelling data

------------------------------------------------------------------------

### Method `create()`

create all directories

#### Usage

    ProjectStructure$create()

------------------------------------------------------------------------

### Method `reset()`

empty modelling_path and qc_path folder.

#### Usage

    ProjectStructure$reset()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ProjectStructure$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
tmp <- ProjectStructure$new("./test_project",
project_Id  = 3000,
order_Id = 6200,
workunit_Id = 23000,
inputAnnotation = ".",
inputData = "."
)
tmp$qc_path()
#> [1] "./test_project/qc_results"
tmp$modelling_path()
#> [1] "./test_project/modelling_results"

tmp$modelling_path()
#> [1] "./test_project/modelling_results"
tmp$modelling_dir
#> [1] "modelling_results"
tmp$modelling_path("second_model")
#> [1] "./test_project/modelling_results" "./test_project/second_model"     
tmp$create()
#> NULL
tmp$reset()
#> NULL
unlink("./test_project", recursive = TRUE)
```
