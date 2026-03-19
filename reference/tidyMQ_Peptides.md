# parse MQ peptides.txt

parse MQ peptides.txt

parse MQ peptides.txt

## Usage

``` r
tidyMQ_Peptides(MQPeptides, proteotypic_only = TRUE)

tidyMQ_Peptides(MQPeptides, proteotypic_only = TRUE)
```

## Arguments

- MQPeptides:

  data.frame generated with read.csv("peptide.txt",sep = "\t",
  stringsAsFactors = FALSE)

## See also

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Evidence()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Evidence.md),
[`tidyMQ_ProteinGroups()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_ProteinGroups.md)

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Evidence()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Evidence.md),
[`tidyMQ_ProteinGroups()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_ProteinGroups.md)

## Examples

``` r

peptide_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")

peptides_txt <- read.csv(
  unz(peptide_txt, "peptides.txt"),
  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mq_peptides <- tidyMQ_Peptides(peptides_txt)



peptide_txt <- prolfqua::find_package_file(
  "prolfquapp", "samples/maxquant_txt/tiny2.zip")

peptides_txt <- read.csv(
  unz(peptide_txt, "peptides.txt"),
  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mq_peptides <- tidyMQ_Peptides(peptides_txt)
```
