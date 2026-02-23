# get_annot_from_fasta

get_annot_from_fasta

## Usage

``` r
get_annot_from_fasta(
  fasta.files,
  pattern_decoys = "^REV_|^rev_",
  isUniprot = TRUE,
  min_length = 7,
  max_length = 30,
  include_seq = FALSE
)
```

## Arguments

- fasta.files:

  path to fasta file(s) or connection

- pattern_decoys:

  regex for decoy sequence IDs

- isUniprot:

  if TRUE parse UniProt-style headers

- min_length:

  minimum tryptic peptide length

- max_length:

  maximum tryptic peptide length

- include_seq:

  if TRUE include protein sequences

## Examples

``` r
fasta_conn <- textConnection(prolfquapp:::.getSequences())
#testthat::expect_error(prolfquapp::get_annot_from_fasta(fasta_conn, pattern_decoys = "" ))
close(fasta_conn)
fasta_conn <- textConnection(prolfquapp:::.getSequences())
prolfquapp::get_annot_from_fasta(fasta_conn, pattern_decoys = "^REV_|^rev" )
#> INFO [2026-02-23 21:06:55] get_annot : finished reading
#> INFO [2026-02-23 21:06:55] get_annot : extract headers
#> INFO [2026-02-23 21:06:55] get_annot : all seq : 11
#> INFO [2026-02-23 21:06:55] removing decoy sequences usin patter : ^REV_|^rev
#> INFO [2026-02-23 21:06:55] get_annot nr seq after decoy removal: 7
#> INFO [2026-02-23 21:06:55] get_annot : isUniprot : TRUE
#> INFO [2026-02-23 21:06:55] get_annot : extracted gene names
#> INFO [2026-02-23 21:06:55] get_annot : protein length
#> INFO [2026-02-23 21:06:55] get_annot : nr of tryptic peptides per protein computed.
#>                                          fasta.id
#> sp|A0A385XJL2|YGDT_ECOLI sp|A0A385XJL2|YGDT_ECOLI
#> sp|A5A615|YNCL_ECOLI         sp|A5A615|YNCL_ECOLI
#> sp|P03018|UVRD_ECOLI         sp|P03018|UVRD_ECOLI
#> sp|P04982|RBSD_ECOLI         sp|P04982|RBSD_ECOLI
#> sp|P04994|EX7L_ECOLI         sp|P04994|EX7L_ECOLI
#> zz|Y-FGCZCont00001|           zz|Y-FGCZCont00001|
#> zz|Y-FGCZCont00002|           zz|Y-FGCZCont00002|
#>                                                                                                              fasta.header
#> sp|A0A385XJL2|YGDT_ECOLI                         Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1
#> sp|A5A615|YNCL_ECOLI             Uncharacterized protein YncL OS=Escherichia coli (strain K12) OX=83333 GN=yncL PE=1 SV=1
#> sp|P03018|UVRD_ECOLI                          DNA helicase II OS=Escherichia coli (strain K12) OX=83333 GN=uvrD PE=1 SV=1
#> sp|P04982|RBSD_ECOLI                        D-ribose pyranase OS=Escherichia coli (strain K12) OX=83333 GN=rbsD PE=1 SV=3
#> sp|P04994|EX7L_ECOLI     Exodeoxyribonuclease 7 large subunit OS=Escherichia coli (strain K12) OX=83333 GN=xseA PE=1 SV=2
#> zz|Y-FGCZCont00001|                                            zz_FGCZCont0000_P61626_LYSC_HUMAN blastpHomologue_5.0e-107
#> zz|Y-FGCZCont00002|                                                 zz_FGCZCont0001_P02534_K1M1_SHEEP blastpHomologue_0.0
#>                              proteinname gene_name protein_length
#> sp|A0A385XJL2|YGDT_ECOLI      A0A385XJL2      ygdT             48
#> sp|A5A615|YNCL_ECOLI              A5A615      yncL             31
#> sp|P03018|UVRD_ECOLI              P03018      uvrD             60
#> sp|P04982|RBSD_ECOLI              P04982      rbsD             60
#> sp|P04994|EX7L_ECOLI              P04994      xseA             60
#> zz|Y-FGCZCont00001|      Y-FGCZCont00001                       60
#> zz|Y-FGCZCont00002|      Y-FGCZCont00002                       60
#>                          nr_tryptic_peptides
#> sp|A0A385XJL2|YGDT_ECOLI                   0
#> sp|A5A615|YNCL_ECOLI                       1
#> sp|P03018|UVRD_ECOLI                       4
#> sp|P04982|RBSD_ECOLI                       3
#> sp|P04994|EX7L_ECOLI                       2
#> zz|Y-FGCZCont00001|                        4
#> zz|Y-FGCZCont00002|                        1
close(fasta_conn)
```
