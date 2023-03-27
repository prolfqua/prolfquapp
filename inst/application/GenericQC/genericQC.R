#R
# 20211210 WEW/CP make it work for WU272669

##### QCs


protein <- read.csv("combined_protein.tsv",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE ) |>
  tidyr::as_tibble() |>
  prolfqua::tidy_FragPipe_combined_protein() |>
  dplyr::filter(combined.total.peptides > 1)


annotation <- try(read.delim("samples.txt"))
if (!inherits(annotation, 'try-error')) {
  annotation$inputresource.name <- tools::file_path_sans_ext(tools::file_path_sans_ext(annotation$inputresource.name))
  annotation$sample.name <- make.unique(annotation$sample.name)

  # unique.spectral.count
  # dataset; adapt name 'group' from dataset
  protein <- dplyr::inner_join(annotation, protein,
                               by = c("inputresource.name" = "raw.file"),
                               multiple = "all")

  # MSFragger specific (moving target)
  atable <- prolfqua::AnalysisTableAnnotation$new()
  atable$fileName = "inputresource.name"
  atable$sampleName = "sample.name"
  atable$hierarchy[["protein_Id"]] <- c("protein")
  atable$hierarchyDepth <- 1
  atable$setWorkIntensity("razor.intensity")
  atable$factors[["group"]] = "groupingvar.name"
  atable$factorDepth <- 1

} else{
  annotation <- data.frame(raw.file = unique(protein$raw.file), group = 'file')
  protein <- dplyr::inner_join(annotation, protein,
                               multiple = "all")

  # MSFragger specific (moving target)
  atable <- prolfqua::AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  atable$hierarchy[["protein_Id"]] <- c("protein")
  atable$hierarchyDepth <- 1
  atable$setWorkIntensity("razor.intensity")
  atable$factors[["group"]] = "group"
  atable$factorDepth <- 1
}

config <- prolfqua::AnalysisConfiguration$new(atable)
adata <- prolfqua::setup_analysis(protein, config)
lfqdata <- prolfqua::LFQData$new(adata, config)
lfqdata$remove_small_intensities(threshold = 1)

ps <- prolfqua::ProjectStructure$new(outpath = ".",
                                     project_Id = "",
                                     workunit_Id = basename(getwd()),
                                     order_Id = "",
                                     inputAnnotation = NULL,
                                     inputData = NULL)

ps$create()

#adapt project/order ID
if(FALSE){
  prolfqua::render_MQSummary_rmd(lfqdata$data,
                                 config$clone(deep = TRUE),
                                 ps, format = "html")
}else {
  prr = list(
    data = lfqdata$data,
    configuration = lfqdata$config,
    project_conf = ps,
    pep = FALSE
  )
  rmarkdown::render("../QCandSSE.Rmd",params = prr)

}
