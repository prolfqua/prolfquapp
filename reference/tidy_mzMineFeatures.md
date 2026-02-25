# convert mzmine features to tidy table

convert mzmine features to tidy table

## Usage

``` r
tidy_mzMineFeatures(data)
```

## Arguments

- data:

  path to csv or data frame of mzMine features

## Examples

``` r
# example code
if(FALSE){
file_path = "outputs-20250407T1707/mzmine/result_features.csv"
raw_df <- readr::read_csv(file_path)
res <- tidy_mzMineFeatures(raw_df)
 file_path <- paste0(
   "/Users/witoldwolski/__checkout/prolfquapp/",
   "inst/application/mzMine/out_results_zip/",
   "mzmine/result_features.csv")
#file_path = "WU323671_mzMine_o35537_WpH9V2_neg_v2_result/mzmine/result_features.csv"
x <- readr::read_csv(file_path)
raw_df
res <- tidy_mzMineFeatures(raw_df)
head(res)

}
```
