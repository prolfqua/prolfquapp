# extract intensities and annotations from MQ proteinGroups.txt

extract intensities and annotations from MQ proteinGroups.txt

extract intensities and annotations from MQ proteinGroups.txt

## Usage

``` r
tidyMQ_ProteinGroups(MQProteinGroups)

tidyMQ_ProteinGroups(MQProteinGroups)
```

## Arguments

- MQProteinGroups:

  data.frame generated with read.csv("peptide.txt",sep="\t",
  stringsAsFactors=FALSE)

## See also

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Evidence()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Evidence.md),
[`tidyMQ_Peptides()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Peptides.md)

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Evidence()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Evidence.md),
[`tidyMQ_Peptides()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Peptides.md)

## Examples

``` r
protein_txt <- prolfqua::find_package_file("prolfquapp","samples/maxquant_txt/tiny2.zip")
protein_txt <- read.csv(
  unz(protein_txt, "proteinGroups.txt"),
  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mq_proteins <-tidyMQ_ProteinGroups(protein_txt)

protein_txt <- prolfqua::find_package_file(
  "prolfquapp", "samples/maxquant_txt/tiny2.zip")
protein_txt <- read.csv(
  unz(protein_txt, "proteinGroups.txt"),
  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mq_proteins <-tidyMQ_ProteinGroups(protein_txt)
```
