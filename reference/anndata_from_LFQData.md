# convert lfqdata to anndata

convert lfqdata to anndata

## Usage

``` r
anndata_from_LFQData(lfqdata, pannot)
```

## Arguments

- lfqdata:

  LFQData object

- pannot:

  ProteinAnnotation object

## Examples

``` r
# example code
library(prolfqua)
library(prolfquapp)
lfqdata <- prolfqua::sim_lfq_data_2Factor_config()
#> creating sampleName from fileName column
#> completing cases
#> completing cases done
#> setup done
lfqdata <- LFQData$new(lfqdata$data, lfqdata$config)
lfqdata$data$protein_Id <- add_RevCon(lfqdata$data$protein_Id)
pids <- grep("^zz|^REV", unique(lfqdata$data$protein_Id), value = TRUE, invert = TRUE)
addannot <- data.frame(
  protein_Id = pids,
  description = stringi::stri_rand_strings(length(pids), 13)
)

addannot <- addannot |> tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
pannot <- ProteinAnnotation$new(lfqdata,
  addannot,
  description = "description",
  cleaned_ids = "cleanID",
  pattern_contaminants = "^zz",
  pattern_decoys = "^REV"
)
#> Warning: no exp_nr_children column specified, computing using nr_obs_experiment function

#debug(anndata_from_LFQData)
anndata_from_LFQData(lfqdata, pannot)
#> converting to layers: abundance, qValue, nr_peptides
#> InMemoryAnnData object with n_obs × n_vars = 16 × 10
#>     obs: 'sample', 'sampleName', 'Treatment', 'Background'
#>     var: 'protein_Id', 'cleanID', 'description', 'nr_peptides'
#>     layers: 'abundance', 'qValue', 'nr_peptides'
#anndataR::write_h5ad(adata, path = "test.h5ad", mode = "w")
```
