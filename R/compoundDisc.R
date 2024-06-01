#' massage CD output compound table.
#' @export
massage_CD <- function(in_file){
  xd <- readxl::read_excel(in_file)
  xd$my_C_ID <- 1:nrow(xd)
  annot <- xd |> dplyr::select("my_C_ID", "Checked","Tags",
                        "Structure","Name","Formula","Annot. Source: Predicted Compositions","Annot. Source: mzCloud Search",
                        "Annot. Source: mzVault Search","Annot. Source: ChemSpider Search","Annot. Source: MassList Search",
                        "Annotation MW", "Calc. MW","m/z","RT [min]")
  colnames(annot) <- gsub("[[:space:].:/]+", "_",colnames(annot))
  colnames(annot) <- gsub("\\[|\\]","",colnames(annot))
  annot <- annot |> dplyr::mutate(FormulaB = stringr::str_replace_all(Formula, " ",""))
  annot <- annot |> tidyr::unite("NewID", c("FormulaB", "m_z", "RT_min"), remove = FALSE)


  tolong <- xd |> dplyr::select("my_C_ID",tidyselect::starts_with(c("Area:","Gap Status:","Gap Fill Status:","Peak Rating:")),)
  xdl <- tolong |> tidyr::pivot_longer(
    cols = tidyselect::starts_with(c("Area:","Gap Status:","Gap Fill Status:","Peak Rating:")),
    names_to = c(".value","filename","file_id"),names_pattern = "(.*)\\: (.*)(\\s\\(F\\d+\\))" )



  colnames(xdl) <- gsub("# ", "", colnames(xdl) )
  colnames(xdl) <- gsub("[[:space:]]","_",colnames(xdl))


  xdl <- xdl |> dplyr::mutate(file_id = gsub(" ","",gsub("\\(|\\)","", file_id)))

  xdl$SampleName <- xdl$file_id

  xdl <- dplyr::inner_join(annot, xdl, by = "my_C_ID")
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
  nr <- sum(annot$File.Name %in% sort(unique(xdl$filename)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(xdl$raw.file)))
  stopifnot( nr > 0)
  atable <- annotation$atable

  atable$ident_Score = "score"
  atable$ident_qValue = "qValue"
  atable$fileName = "File.Name"
  atable$sampleName = "SampleName"
  atable$hierarchy[["metabolite_Id"]] <- c("NewID")
  atable$set_response("Area")

  byv <- c("filename")

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
