#' massage CD output compound table.
#' @export
massage_CD <- function(in_file ){

  xd <- if (is.character(in_file) && file.exists(in_file)) {
    readxl::read_excel(in_file)
  } else if (is.data.frame(in_file)) {
    in_file
  } else { stopifnot("expecting data frame or path got : ", class(in_file))}
  xd$my_C_ID <- 1:nrow(xd)
  annot <- xd |> dplyr::select(any_of(c("my_C_ID", "Checked","Tags",
                                        "Structure","description" = "Name","Formula","Annot Source Predicted Compositions","Annot Source mzCloud Search",
                                        "Annot Source mzVault Search","Annot Source ChemSpider Search","Annot Source MassList Search",
                                        "Annotation MW", "Calc MW","mz",RT_min = "RT in min")))
  colnames(annot) <- gsub("[[:space:].:/]+", "_",colnames(annot))
  colnames(annot) <- gsub("\\[|\\]","",colnames(annot))
  annot <- annot |> dplyr::mutate(FormulaB = stringr::str_replace_all(Formula, " ",""))
  annot <- annot |> tidyr::unite(
    "metabolite_feature_Id",
    c("my_C_ID","FormulaB", "mz", "RT_min"), sep = "_", remove = FALSE)

  columns <- c("Area","Gap Status","Gap Fill Status","Peak Rating")

  tolong <- xd |> dplyr::select("my_C_ID",tidyselect::starts_with(columns))
  tolong <- select(tolong, -all_of(c("Area Max", "Area SD", "Area CV in Percent")))

  sum(grepl(columns[1], colnames(tolong)))
  sum(grepl(columns[2], colnames(tolong)))
  sum(grepl(columns[3], colnames(tolong)))
  sum(grepl(columns[4], colnames(tolong)))

  xdl <- tolong |> tidyr::pivot_longer(
    cols = tidyselect::starts_with(c("Area","Gap Status","Gap Fill Status","Peak Rating")),
    names_to = c(".value","filename","file_id"),names_pattern = "(.*)\\s(.*)\\s(F\\d+)" )
  xdl$`Gap Fill Status` |> table()

  colnames(xdl) <- gsub("# ", "", colnames(xdl) )
  colnames(xdl) <- gsub("[[:space:]]","_",colnames(xdl))
  xdl <- xdl |> dplyr::mutate(file_id = gsub(" ","",gsub("\\(|\\)","", file_id)))

  # use nr_children to encode gap status.
  xdl <- xdl |> dplyr::mutate(
    nr_children = dplyr::case_when(
      Gap_Status == "Full gap" ~ 0,
      Gap_Status == "Missing ions" ~ 1,
      Gap_Status == "No gap" ~ 2,
      TRUE ~ 3))

  xdl <- dplyr::inner_join(annot, xdl, by = "my_C_ID")
  return(xdl)
}

#' load compound discoverer (CD) files
#' @param in_file excel file produced by CD
#' @param annotation list returned by `read_annotation` function
#' @export
preprocess_CD <- function(
    xdl,
    annotation,
    .func_massage = prolfquapp::massage_CD
){
  annot <- annotation$annot
  nr <- sum(annot$File.Name %in% sort(unique(xdl$filename)))
  logger::log_info("nr : ", nr, " files annotated out of ", length(unique(xdl$filename)))
  stopifnot( nr > 0)


  atable <- annotation$atable
  atable$sampleName = "file_id"
  atable$hierarchy[["metabolite_feature_Id"]] <- c("metabolite_feature_Id")
  atable$set_response("Area")
  byv <- c("filename")
  names(byv) <- atable$fileName
  byv <- c(byv, intersect(colnames(annot), colnames(xdl)))

  peptide <- dplyr::inner_join(annot, xdl, by = byv, multiple = "all")
  config <- prolfqua::AnalysisConfiguration$new(atable)
  adata <- prolfqua::setup_analysis(peptide, config)
  lfqdata <- prolfqua::LFQData$new(adata, config)
  lfqdata$remove_small_intensities()

  m_annot <- xdl |>
    dplyr::select("metabolite_feature_Id", "Checked", "Tags", "Structure", "description","Formula", starts_with("Annot_")) |>
    dplyr::distinct()
  # handle not identified
  m_annot$nr_compounds <- ifelse(m_annot$Checked, 2 ,1)

  m_annot <- m_annot |> dplyr::mutate(IDcolumn = metabolite_feature_Id)
  prot_annot <- prolfquapp::ProteinAnnotation$new(
    lfqdata , m_annot, description = "description",
    cleaned_ids = "IDcolumn",
    full_id = "Checked",
    exp_nr_children = "nr_compounds",
    pattern_contaminants = "FALSE",
    pattern_decoys = NULL
  )
  return(list(lfqdata = lfqdata , protein_annotation = prot_annot))
}



