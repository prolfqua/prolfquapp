# Write dataset to file in csv, tsv, or xlsx format

Write dataset to file in csv, tsv, or xlsx format

## Usage

``` r
write_annotation_file(data, file_path)
```

## Arguments

- data:

  data frame to write

- file_path:

  output file path (csv, tsv, or xlsx)

## Examples

``` r
ds <- data.frame(channel = c("A","B","C"), Name = NA, Subject = NA, Group = NA, Control = NA)
write_annotation_file(ds, file_path = file.path(tempdir(),"test.xlsx"))
```
