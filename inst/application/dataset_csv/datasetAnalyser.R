dataset_csv_examples <- system.file("application/dataset_csv",package = "prolfquapp")
files <- dir(dataset_csv_examples, pattern = "*.csv")

res <- list()
for (i in seq_along(files)) {
  print(i)
  res[[i]] <- readr::read_csv(files[i])
}

for (i in seq_along(res)) {
  atable <- prolfqua::AnalysisConfiguration$new()
  atable$fileName = "channel"
  atable$hierarchy[["protein_Id"]] <- c("Protein")
  atable$hierarchy[["peptide_Id"]] <- c("Peptide")
  tmp <- prolfquapp::dataset_set_factors(atable, res[[i]] )
  cat(i, " : " , length(tmp$atable$factors), "factors : ", paste(tmp$atable$factors, collapse = "; "), "\n")
}
