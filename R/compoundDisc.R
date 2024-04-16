#' massage CD output compound table.
#' @export
massage_CD <- function(in_file){
  xd <- readxl::read_excel(in_file)
  # xd <- xd |> dplyr::filter(!is.na(Name), Name != "Tags", Formula != "Checked", Formula != FALSE)
  xd <- xd |> dplyr::select(1:max(grep("^Area",colnames(xd))))
  grep("Area",colnames(xd), value = TRUE)
  xd$`Area (Max.)` <- NULL
  xd <- xd |> dplyr::select(-starts_with("Group"),
                            -starts_with("Ratio"),
                            -starts_with("Log2"),
                            -starts_with("P-value"))

  xd <- xd |> dplyr::mutate(FormulaB = stringr::str_replace_all(Formula, " ",""))
  xd <- xd |> tidyr::unite("NewID", c("FormulaB", "m/z", "RT [min]"), remove = FALSE)
  xdl <- xd |> tidyr::pivot_longer(cols = tidyselect::starts_with("Area"),names_to = "Sample", values_to = "Area")
  colnames(xdl) <- gsub("# ", "", colnames(xdl) )
  colnames(xdl) <- gsub("Annot. Source:", "Annot", colnames(xdl))
  colnames(xdl) <- make.names(colnames(xdl))


  xdl <- xdl |> dplyr::mutate(SampleS = gsub("Area: ","", Sample))
  xdl <- xdl |> dplyr::mutate(SampleS = gsub("\\(|\\)","", SampleS))
  xdl <- xdl |> tidyr::separate(SampleS, c("raw.file", "StudyFile"), sep=" ")
  xdl$StudyFile |> unique()
  xdl$compound_name <- xdl$Name
  xdl$SampleName <- xdl$StudyFile
  xdl <- xdl
  xdl$score <- 10
  xdl$qValue <- 0
  xdl$nr_compounds <- 1
  return(xdl)
}

#' load compound discoverer (CD) files
#' @param in_file excel file produced by CD
#' @param annotation list returned by `read_annotation` function
#' @export
preprocess_CD <- function(in_file,
                          annotation,
                          .func_massage = prolfquapp::massage_CD
){

  xdl <- .func_massage(in_file)


  annot <- annotation$annot
  nr <- sum(annot$File.Name %in% sort(unique(xdl$raw.file)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(xdl$raw.file)))
  stopifnot( nr > 0)
  atable <- annotation$atable
  atable$ident_Score = "score"
  atable$ident_qValue = "qValue"
  atable$fileName = "File.Name"
  atable$sampleName = "SampleName"
  atable$hierarchy[["metabolite_Id"]] <- c("NewID")
  atable$set_response("Area")

  byv <- c("raw.file")

  names(byv) <- atable$fileName
  peptide <- dplyr::inner_join(annot, xdl, by = byv, multiple = "all")

  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  prot_annot <- prolfquapp::build_protein_annot(
    lfqdata,
    xdl,idcol = c("metabolite_Id" = "NewID"),
    cleaned_protein_id = "NewID",
    protein_description = "compound_name",
    nr_children = "nr_compounds",
    more_columns = c("Formula","mzVault.Results", "ChemSpider.Results", "mzCloud.Results" )
  )

  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}
