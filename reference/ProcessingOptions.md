# processing options R6 class

processing options R6 class

processing options R6 class

## See also

Other ProlfquAppConfig:
[`ExternalReader`](https://prolfqua.github.io/prolfquapp/reference/ExternalReader.md),
[`ProjectSpec`](https://prolfqua.github.io/prolfquapp/reference/ProjectSpec.md),
[`ProlfquAppConfig`](https://prolfqua.github.io/prolfquapp/reference/ProlfquAppConfig.md),
[`make_DEA_config_R6()`](https://prolfqua.github.io/prolfquapp/reference/make_DEA_config_R6.md),
[`set_list_to_R6()`](https://prolfqua.github.io/prolfquapp/reference/set_list_to_R6.md)

## Public fields

- `transform`:

  data transformation method

- `aggregate`:

  protein abundance estimation method

- `diff_threshold`:

  difference threshold

- `FDR_threshold`:

  FDR threshold

- `remove_cont`:

  should contaminants be removed

- `remove_decoys`:

  should decoys be removed

- `pattern_decoys`:

  decoy patterns

- `pattern_contaminants`:

  pattern contaminants

- `nr_peptides`:

  number of peptides

- `interaction`:

  model with interactions default FALSE

- `model_missing`:

  model missigness, default TRUE

- `model`:

  name of the model to use "prolfqua", "SE", "ROPECA", default
  "prolfqua"

- `other`:

  list with additional options

- `internal`:

  list of internal calibrants

## Methods

### Public methods

- [`ProcessingOptions$clone()`](#method-ProcessingOptions-clone)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ProcessingOptions$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
