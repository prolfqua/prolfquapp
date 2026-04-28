# Decorates LFQData with a row annotation and some protein specific functions.

Decorates LFQData with a row annotation and some protein specific
functions.

Decorates LFQData with a row annotation and some protein specific
functions.

## Public fields

- `row_annot`:

  data.frame containing further information

- `pID`:

  column with protein ids

- `full_id`:

  column with protein id e.g. sp\| can be same as pID

- `description`:

  name of column containing descriptions

- `cleaned_ids`:

  vector with columns containing addition IDs

- `exp_nr_children`:

  name of columns with the number of peptides

- `pattern_contaminants`:

  pattern_contaminants

- `pattern_decoys`:

  pattern_decoys

## Methods

### Public methods

- [`ProteinAnnotation$new()`](#method-ProteinAnnotation-new)

- [`ProteinAnnotation$annotate_decoys()`](#method-ProteinAnnotation-annotate_decoys)

- [`ProteinAnnotation$annotate_contaminants()`](#method-ProteinAnnotation-annotate_contaminants)

- [`ProteinAnnotation$get_summary()`](#method-ProteinAnnotation-get_summary)

- [`ProteinAnnotation$nr_clean()`](#method-ProteinAnnotation-nr_clean)

- [`ProteinAnnotation$clean()`](#method-ProteinAnnotation-clean)

- [`ProteinAnnotation$filter_by_nr_children()`](#method-ProteinAnnotation-filter_by_nr_children)

- [`ProteinAnnotation$clone()`](#method-ProteinAnnotation-clone)

------------------------------------------------------------------------

### Method `new()`

initialize

#### Usage

    ProteinAnnotation$new(
      lfqdata,
      row_annot = NULL,
      description = NULL,
      cleaned_ids = NULL,
      full_id = NULL,
      exp_nr_children = "nr_peptides",
      pattern_contaminants = NULL,
      pattern_decoys = NULL
    )

#### Arguments

- `lfqdata`:

  data frame from
  [`setup_analysis`](https://wolski.github.io/prolfqua/reference/setup_analysis.html)

- `row_annot`:

  data frame with row annotation. Must have columns matching
  `config$hierarchy_keys_depth()`

- `description`:

  name of column with description

- `cleaned_ids`:

  names of columns with cleaned Ids

- `full_id`:

  column with full protein ID

- `exp_nr_children`:

  column with the number of children

- `pattern_contaminants`:

  pattern_contaminants

- `pattern_decoys`:

  pattern_decoys

------------------------------------------------------------------------

### Method `annotate_decoys()`

annotate rev sequences

#### Usage

    ProteinAnnotation$annotate_decoys()

#### Arguments

- `pattern`:

  default "REV\_"

------------------------------------------------------------------------

### Method `annotate_contaminants()`

annotate contaminants

#### Usage

    ProteinAnnotation$annotate_contaminants()

#### Arguments

- `pattern`:

  default "^zz\|^CON"

------------------------------------------------------------------------

### Method `get_summary()`

get summary

#### Usage

    ProteinAnnotation$get_summary()

------------------------------------------------------------------------

### Method `nr_clean()`

get number of neither contaminants nor decoys

#### Usage

    ProteinAnnotation$nr_clean(contaminants = TRUE, decoys = TRUE)

#### Arguments

- `contaminants`:

  remove contaminants

- `decoys`:

  remove decoys return number of cleans

------------------------------------------------------------------------

### Method `clean()`

remove REV and CON sequences

#### Usage

    ProteinAnnotation$clean(contaminants = TRUE, decoys = TRUE)

#### Arguments

- `contaminants`:

  remove contaminants

- `decoys`:

  remove decoys

------------------------------------------------------------------------

### Method `filter_by_nr_children()`

filter by number children

#### Usage

    ProteinAnnotation$filter_by_nr_children(exp_nr_children = 2)

#### Arguments

- `exp_nr_children`:

  minimum number of children required

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ProteinAnnotation$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#> creating sampleName from file_name column
#> completing cases
#> completing cases done
#> setup done
lfq0 <- prolfqua::LFQData$new(istar$data, istar$config)
xd1 <- prolfqua::nr_children_experiment(lfq0$data_long(), lfq0$response(),
  lfq0$relevant_hierarchy_keys(), lfq0$file_name(), lfq0$nr_children_col())

xd2 <- prolfqua::nr_features_experiment(lfq0$data_long(), lfq0$hierarchy_keys(),
  lfq0$relevant_hierarchy_keys())
xd1$nr_child_exp |> table()
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 
#> 22 26 16 10  9  3  5  4  1  1  1  2 

lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
tmp <- lfqdata$data_long()
tmp$protein_Id <- add_RevCon(tmp$protein_Id)
lfqdata$set_data(tmp)
pids <- grep("^zz|^REV", unique(lfqdata$data_long()$protein_Id), value = TRUE, invert = TRUE)
addannot <- data.frame(
  protein_Id = pids,
  description = stringi::stri_rand_strings(length(pids), 13)
)

addannot <- addannot |> tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
# ProteinAnnotation$debug("initialize")
# debug(nr_obs_sample)
xd4 <- prolfqua::nr_obs_sample(lfqdata$data_long(), lfqdata$response(),
  lfqdata$relevant_hierarchy_keys(), lfqdata$file_name(), lfqdata$nr_children_col())
xd3 <- prolfqua::nr_features_experiment(lfqdata$data_long(), lfqdata$hierarchy_keys(),
  lfqdata$relevant_hierarchy_keys())

pannot <- ProteinAnnotation$new(lfqdata,
  addannot,
  description = "description",
  cleaned_ids = "cleanID",
  pattern_contaminants = "^zz",
  pattern_decoys = "^REV"
)
#> Warning: no exp_nr_children column specified, computing using nr_children_experiment
stopifnot(pannot$annotate_decoys() == 10)
stopifnot(pannot$annotate_contaminants() == 5)
dd <- pannot$clean()
pannot$nr_clean()
#> [1] 85
pannot$get_summary()
#>   totalNrOfProteins percentOfContaminants percentOfFalsePositives
#> 1               100                     5                      10
#>   NrOfProteinsNoDecoys
#> 1                   85
stopifnot(nrow(dd) == 85)
tmp <- lfqdata$get_subset(dd)
#> Joining with `by = join_by(protein_Id)`
dx <- pannot$clean(contaminants = TRUE, decoys = FALSE)
stopifnot(nrow(dx) == 95)
dx <- pannot$clean(contaminants = FALSE, decoys = TRUE)
stopifnot(nrow(dx) == 90)
dx2 <- pannot$filter_by_nr_children(exp_nr_children = 2)
dx3 <- pannot$filter_by_nr_children(exp_nr_children = 3)
stopifnot(nrow(dx2) >= nrow(dx3))
```
