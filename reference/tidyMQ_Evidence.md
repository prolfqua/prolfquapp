# read evidence file

read evidence file

read evidence file

## Usage

``` r
tidyMQ_Evidence(Evidence)

tidyMQ_Evidence(Evidence)
```

## Arguments

- Evidence:

  MQ evidence file or zip archive with evidence file

## See also

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Peptides()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Peptides.md),
[`tidyMQ_ProteinGroups()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_ProteinGroups.md)

Other MaxQuant:
[`MaxQuant`](https://prolfqua.github.io/prolfquapp/reference/MaxQuant.md),
[`tidyMQ_Peptides()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_Peptides.md),
[`tidyMQ_ProteinGroups()`](https://prolfqua.github.io/prolfquapp/reference/tidyMQ_ProteinGroups.md)

## Examples

``` r
evidence_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
evidence_txt <- read.csv(unz(evidence_txt,"evidence.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
mq_evidence <- tidyMQ_Evidence(evidence_txt)

evidence_txt <- prolfqua::find_package_file("prolfquapp", "samples/maxquant_txt/tiny2.zip")
evidence_txt <- read.csv(unz(evidence_txt,"evidence.txt"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
mq_evidence <- tidyMQ_Evidence(evidence_txt)
```
