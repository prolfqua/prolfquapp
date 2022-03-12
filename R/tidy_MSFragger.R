tidy_MSFragger_combined_protein_V16 <- function(
  combprot,
  as_list = FALSE
) {
  spcnames = c("Total Spectral Count",
               "Unique Spectral Count",
               "Razor Spectral Count")
  intnames = c("Total Intensity",
               "Unique Intensity",
               "Razor Intensity")
  maxlfqnames = c("Maxlfq Total Intensity",
                  "Maxlfq Unique Intensity",
                  "Maxlfq Razor Intensity")

  protIDcol = "Protein"

  if (is.character(combprot) && file.exists(combprot)) {
    Cprotein <- as_tibble(read.csv(combprot,
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))

  } else if ("tbl_df" %in% class(combprot)) {
    Cprotein <- combprot
  } else {
    stop(class(combprot), " not supported.")
  }

  cnam <- gsub("Total Razor ", "Total ",
               gsub("Unique Razor ","Unique ",
                    gsub(" Intensity$"," Razor Intensity",
                         gsub(" Spectral Count$"," Razor Spectral Count",colnames(Cprotein))
                    )
               )
  )

  ### start processing
  cnam <- cnam
  colnames(Cprotein) <- cnam
  cnam <- cnam[1:which(cnam == "Combined Total Spectral Count")]

  message("annotation columns : ", paste(cnam, collapse = "\n"))
  annot <- Cprotein |> dplyr::select(all_of(cnam))
  colnames(Cprotein)

  extractDataLong <- function(Cprotein, what = "Total Intensity", butNot = NULL){
    cols <- colnames(Cprotein)
    cols <- setdiff( grep(paste0(what,"$"), cols, value = TRUE) , if (is.null(butNot)) {NULL} else { grep(butNot, cols, value = TRUE) })
    gg <- Cprotein |> dplyr::select( all_of(protIDcol), all_of(cols) )

    gg <- gg |> tidyr::pivot_longer(cols = dplyr::ends_with(what), names_to = "raw.file",values_to = what)
    gg <- gg |> dplyr::mutate(raw.file = gsub(paste0(".",what,"$"),"", .data$raw.file))
    gg
  }

  res <- vector( mode = "list", length = length(c(intnames, spcnames)))
  names(res)  <- c(intnames, spcnames)

  for (i in seq_along(c(intnames, spcnames))) {
    message("DD: ", c(intnames, spcnames)[i] )
    res[[c(intnames, spcnames)[i]]] <- extractDataLong(Cprotein, what = c(intnames, spcnames)[i], butNot = "maxlfq" )
  }

  if (sum(grepl(".maxlfq.", colnames(Cprotein))) > 0) {
    res_maxlfq <- vector( mode = "list", length(maxlfqnames))
    names(res_maxlfq)  <- maxlfqnames
    for (i in seq_along(maxlfqnames)) {
      message("DD: ", maxlfqnames[i] )
      res_maxlfq[[maxlfqnames[i] ]] <-  extractDataLong(Cprotein, what = maxlfqnames[i], butNot = NULL )
    }
    res <- c(res, res_maxlfq)
  }

  if (as_list) {
    return( res )
  }

  merged <- Reduce( inner_join , res )
  merged <- inner_join( annot, merged )
  colnames(merged) <- tolower(make.names(colnames(merged)))
  return( merged )
}

