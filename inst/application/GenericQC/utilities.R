find_default_files <- function(path){
  mdir <- function(path, pattern){
    file_info <- file.info(path)
    if (file_info$isdir[1]) {
      res <- dir(path, pattern, full.names = TRUE, recursive = TRUE)
    } else {
      stop("unsupported path. \n")
    }
    return(res)
  }

  fasta.files <- mdir(path, "*.fasta$")
  if (any(grepl("database[0-9]*.fasta$", fasta.files))) {
    fasta.files <- grep("database[0-9]*.fasta$", fasta.files, value = TRUE)
  }
  logger::log_info(paste(fasta.files, collapse = "; "))
  diann.output <- mdir(path,
                       pattern = "report\\.tsv$|diann-output\\.tsv$")
  logger::log_info(diann.output)
  dataset.csv <- mdir(path,
                      pattern = "dataset.csv$")
  res <- list(dataset = dataset.csv, data = diann.output, fasta = fasta.files)
}
