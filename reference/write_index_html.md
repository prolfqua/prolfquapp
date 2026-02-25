# write index.html file with links to all relevant files:

write index.html file with links to all relevant files:

## Usage

``` r
write_index_html(file_path_list, result_dir)
```

## Arguments

- file_path_list:

  named list of output file paths

- result_dir:

  directory for the index.html output

## Examples

``` r
.resdir <- "."
write_index_html(prolfquapp:::.test_links,tempdir())
#> Wrote HTML index to: /tmp/RtmpPMLQDa/index.html
```
