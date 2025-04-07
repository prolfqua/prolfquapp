library(tidyr)
library(dplyr)


tidy_mzMineFeatures <- function(file_path) {
  x <- readr::read_csv(file_path)
  xdrop <- x |> dplyr::select(starts_with("alignment_scores"),
                              starts_with("ion_identities"), starts_with("compound_db_identity"), starts_with("lipid_annotations"), starts_with("molecular_networking"))
  stopifnot(nrow(na.omit(xdrop)) == 0)

  x <- x |> dplyr::select(-starts_with("alignment_scores"),-starts_with("ion_identities"), -starts_with("compound_db_identity"), -starts_with("lipid_annotations"), -starts_with("molecular_networking"))

  colnames(x) <- gsub("^datafile:", "datafile_", colnames(x))
  paste0("feature_",colnames(x)[!grepl("^datafile_",colnames(x))])
  colnames(x)[!grepl("^datafile_",colnames(x))] <- paste0("feature_",colnames(x)[!grepl("^datafile_",colnames(x))])

  xl <- x |> pivot_longer(cols = starts_with("datafile_"),
                          names_to = c("datafile", ".value"),
                          names_sep = "\\.raw:" )
  colnames(xl) <- gsub(":","_", colnames(xl))

  return(xl)
}

# file_path = "outputs-20250407T1707/mzmine/result_features.csv"
# res <- tidy_mzMineFeatures(file_path)

