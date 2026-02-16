# get functions for creating datasets

get functions for creating datasets

## Usage

``` r
dataset_get_functions(preprocess_functions)
```

## Arguments

- preprocess_functions:

  list with get_files, dataset, extra_args entries

## Examples

``` r
# Get dataset functions for DIANN
preprocess_functions <- prolfquapp::prolfqua_preprocess_functions[["DIANN"]]
dataset_funcs <- dataset_get_functions(preprocess_functions)

# Use the functions
if (FALSE) { # \dontrun{
files <- dataset_funcs$files_fn("path/to/data")
dataset <- dataset_funcs$dataset_fn(files, "output_file.csv")
} # }
```
