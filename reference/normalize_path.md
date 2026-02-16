# Function to normalize paths for both Windows and Linux

Function to normalize paths for both Windows and Linux

## Usage

``` r
normalize_path(paths, os = .Platform$OS.type)
```

## Arguments

- paths:

  character vector of file paths

- os:

  operating system type

## Value

normalized path

## Examples

``` r
exp_paths <-c("E:\\projects\\p29033\\TKOiWAT\\20240123_015_S629149_iWAT_FL1.d",
  "E:\\projects\\p29033\\TKOiWAT\\20240123_016_S629150_iWAT_FL2.d",
  "E:\\projects\\p29033\\TKOiWAT\\20240123_017_S629151_iWAT_FL3.d")

normalize_path(exp_paths)
#> [1] "E:/projects/p29033/TKOiWAT/20240123_015_S629149_iWAT_FL1.d"
#> [2] "E:/projects/p29033/TKOiWAT/20240123_016_S629150_iWAT_FL2.d"
#> [3] "E:/projects/p29033/TKOiWAT/20240123_017_S629151_iWAT_FL3.d"
```
