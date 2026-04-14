params <-
list(deanalyse = prolfquapp::example_deanalyse())

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)
dea <- eval(params$deanalyse)
cfg <- dea$prolfq_app_config


## ----samples, eval=TRUE-------------------------------------------------------

tab <- data.frame(table(dea$lfq_data_raw$factors()[dea$lfq_data_raw$relevant_factor_keys()]))
colnames(tab)[length(colnames(tab))] <- "# samples"
knitr::kable(tab, caption = "Nr of samples assigned to each group.")


## ----annotation, eval=TRUE----------------------------------------------------
ld <- dea$lfq_data_raw
fx <- ld$factors()
sr <- ld$get_Summariser()
xx <- sr$hierarchy_counts_sample(nr_children = 1)
xx$isotopeLabel <- NULL
xx <- xx |> dplyr::rename(nr_1 = ld$hierarchy_keys())
x2 <- sr$hierarchy_counts_sample(nr_children = 2)
x2$isotopeLabel <- NULL
x2 <- x2 |> dplyr::rename(nr_2 = ld$hierarchy_keys())

fx <- dplyr::left_join(fx,xx, by = ld$sample_name())
fx <- dplyr::left_join(fx,x2, by = ld$sample_name())


knitr::kable(fx, caption = paste0(
  "LC-MS samples annotation table. The content of the sampleName column ",
  "is used as a short form to plot labels. The group to which a sample ",
  "is assigned is shown in the column group."))


## ----nrPerSample, fig.cap="Number of identified proteins across samples.", fig.with=10, fig.height=7----
sum <- dea$lfq_data_raw$get_Summariser()
sum$plot_hierarchy_counts_sample(nr_children = 1)


## ----nrPerSample2, fig.cap="Number of identified proteins across samples.", fig.with=10, fig.height=7----
sum <- dea$lfq_data_raw$get_Summariser()
sum$plot_hierarchy_counts_sample(nr_children = 2)


## ----countProtWithNAs---------------------------------------------------------
res <- dea$lfq_data_raw$to_wide(as.matrix = TRUE)$data
res[!is.na(res)] <- 0
res[is.na(res)] <- 1
allrows <- nrow(res)
res <- res[apply(res,1, sum) > 0, , drop = FALSE]
allNas <- nrow(res)

## ----prepHeatmap--------------------------------------------------------------
pl <- dea$lfq_data_raw$get_Plotter()
nah <- pl$na_heatmap()



## ----naHeat, fig.width=7, fig.height=6, dpi=300, fig.cap="(ref:naHeat)", fig.alt="", eval=TRUE----
nah

## ----vennProteins, fig.cap="(ref:vennProteins)",fig.alt=""--------------------
pups <- prolfqua::upset_interaction_missing_stats(dea$lfq_data_raw, tr = 1)
UpSetR::upset(pups$data, order.by = "freq", nsets = pups$nsets)


## ----normalization, results = 'asis'------------------------------------------
if (cfg$processing_options$transform == "robscale") {
  cat(
    "To do this the z-score of the $\\log_2$ transformed protein",
    "abundances are computed.",
    "Because we need to estimate the protein differences on the",
    "original scale, we have to multiply the $z$-score by the average",
    "standard deviation of all the $N$ samples in the experiment.",
    "After normalization all samples have an equal mean and variance",
    "and a similar distribution.")
} else if (cfg$processing_options$transform == "vsn") {
  cat("To do this the variance stabilizing normalization (vsn) was applied [@HuberVSN2002].")
} else if (cfg$processing_options$transform == "none") {
  cat(
    "However, in some circumstances it is advisable not to normalize",
    "the data, e.g. in case of affinity purification experiments,",
    "or when the requirements of sufficient similarity among samples",
    "are not met.")
}


## ----normalized,  fig.width=8, fig.height=5,dpi=300, fig.cap="(ref:normalized)",fig.alt="", eval=TRUE----

pl <- dea$lfq_data_raw$get_Plotter()

showlegend <- if (nrow(dea$lfq_data_raw$factors()) > 12) {FALSE} else {TRUE}
p1 <- pl$intensity_distribution_density(legend = showlegend)

plTransformed <- dea$lfq_data$get_Plotter()
p2 <- plTransformed$intensity_distribution_density(legend = showlegend) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

gridExtra::grid.arrange(p1, p2, nrow = 1)


## ----internalCalibration, results = 'asis'------------------------------------
internal_ids <- cfg$processing_options$internal
if (length(internal_ids) > 0 && !all(is.na(internal_ids)) && !all(internal_ids == "")) {
  cat("In addition to the global normalization described above, **internal standard normalization** was applied.")
  cat(" The following protein(s) were used as internal reference standard(s):\n\n")
  cat(paste0("- `", internal_ids, "`", collapse = "\n"), "\n\n")
  cat("After transformation, the per-sample intensity of each reference",
    "protein was subtracted from all proteins in that sample",
    "(`center_to_reference`). ")
  cat("As a result, the reference protein(s) have an abundance of exactly 0 in the normalized data, ")
  cat("and all other protein abundances are expressed relative to the reference.\n\n")
}

## ----SDViolin, fig.height=3, fig.width=7, fig.cap="(ref:SDViolin)",fig.alt="", eval=TRUE----
stR <- dea$lfq_data_raw$get_Stats()
pA <- stR$violin() + ggplot2::theme_classic() + ggplot2::labs(tag = "A")
pA


## ----CVtable, eval=TRUE-------------------------------------------------------
resR <- stR$stats_quantiles()
C <- dplyr::bind_rows(CV = resR$wide  |> dplyr::filter(probs == 0.5) |> round(digits = 2))
C <- C |> dplyr::mutate( what = c("CV"), .before = 1 )
C$probs <- NULL
knitr::kable(C, caption = "Median of coefficient of variation (CV).")


## ----generateHeatmaps, inlcude = FALSE----------------------------------------
pl <- dea$lfq_data$get_Plotter()
ph <- pl$heatmap()


## ----heatmap, fig.width=7, fig.height=8, dpi=300, fig.cap="(ref:heatmap)",fig.alt="", eval=TRUE----
ph

## ----pca, fig.cap = "(ref:pca)" ,fig.alt="", fig.width=7, fig.height=7--------
dea$lfq_data$get_Plotter()$pca_plotly()

## ----contrtable---------------------------------------------------------------
df <- data.frame(name = names(dea$contrasts), contrast = dea$contrasts)
rownames(df) <- NULL
knitr::kable(df, caption = "Name of difference (Contrast), and formula used to compute it.")


## -----------------------------------------------------------------------------
hkey <- dea$lfq_data$hierarchy_keys()

## ----setupVolcano, echo = FALSE-----------------------------------------------
datax <- as.data.frame(dea$annotated_contrasts)

## ----tableAllProt, fig.cap="Differential expression analysis results of all proteins."----

datax <- datax |> tidyr::unite("ID" , tidyselect::all_of(hkey), sep = "~")

bb <- datax |> dplyr::select(
  tidyselect::all_of(c("ID",
           "description",
           "contrast",
           "modelName",
           "diff",
           "FDR" )))

bb <- crosstalk::SharedData$new(bb,  as.formula("~ ID"), group = "BB")

DT::datatable(bb, filter = "bottom",
  extensions = "Scroller",
  style = "auto",
  class = "compact",
  options = list(deferRender = TRUE,
                 scrollY = 300,
                 scrollX = 400,
                 scroller = TRUE) ) |>
  DT::formatRound(columns = c("FDR", "diff"), digits = 3)


## ----volcanoplot, fig.cap = "(ref:volcanoplot)",fig.alt="", fig.width=9, fig.height=7, include = TRUE, eval = TRUE----

palette <- c("black","orange")
palette <- setNames(palette, c("Linear_Model_moderated", "Imputed_Mean_moderated"))
xd <- prolfqua::volcano_plotly(
  datax ,
  proteinID = "ID",
  effect = "diff",
  significance = "FDR",
  contrast = "contrast",
  color = "modelName",
  palette = palette,
  xintercept =  c(-cfg$processing_options$diff_threshold,cfg$processing_options$diff_threshold),
  yintercept = cfg$processing_options$FDR_threshold,
  title_size = 10,
  group = "BB")

xd <- lapply(xd, plotly::highlight , off = "plotly_doubleclick")
nrow <- ceiling(length(xd) / 4)
plotly::subplot(xd, shareX = TRUE, shareY = TRUE, nrows = nrow)


## ----nrsignificant, results="markup", eval=TRUE-------------------------------
library(dplyr)
mx <- datax |>
  dplyr::mutate(
    passes = abs(diff) > cfg$processing_options$diff_threshold &
      FDR < cfg$processing_options$FDR_threshold)
x <- mx |>
  dplyr::group_by(contrast) |>
  summarize(
    n = n(),
    Significant = sum(passes, na.rm = TRUE),
    "Not Significant" = n() - sum(passes, na.rm = TRUE))

mycap <- paste0("Number of not significant and significant proteins."  )
knitr::kable(x, caption = mycap)


## ----areThereSig--------------------------------------------------------------
datax_signif <- dea$annotated_contrasts_signif
showSignificant <- TRUE
if (nrow(datax_signif) == 0) {
  showSignificant <- FALSE
}


## ----results = 'asis', eval = showSignificant---------------------------------
cat("The table shown in Figure \\@ref(fig:SigPrey) lists all the significant proteins.")

## ----SigPrey, fig.cap= "(ref:SigPrey)",fig.alt="", eval = showSignificant-----
ctdata <- datax_signif |> tidyr::unite("ID", tidyselect::all_of(hkey), sep = "~")
ctdata <- ctdata |>
  dplyr::select(tidyselect::all_of(c(
    "ID", "description", "contrast",
    "modelName", "diff", "FDR"))) |>
  as.data.frame()

sig <- crosstalk::SharedData$new(ctdata,  as.formula( "~ ID"), group = "BB")
DT::datatable(sig, filter = "bottom",
  extensions = "Scroller",
  style = "auto",
  class = "compact",
  options = list(deferRender = TRUE,
                 scrollY = 300,
                 scrollX = 400,
                 scroller = TRUE) ) |>
  DT::formatRound(columns = c("FDR", "diff"), digits = 3)

## ----prepareHeatmap, eval = showSignificant-----------------------------------
signif <- dea$lfq_data$get_copy()
signif <- signif$get_subset(datax_signif)
if (signif$hierarchy_counts()[2] > 30) {
  rownames = FALSE
} else {
  rownames = TRUE
}
sigheat <- signif$get_Plotter()$raster(rownames = rownames)


## ----makeText, eval = showSignificant, results = 'asis'-----------------------

cat(paste(
  "Furthermore, Figure \\@ref(fig:sigroteins) shows a heatmap",
  "of log2 transformed protein abundances of all",
  "significant calls."))


## ----sigroteins, fig.cap="(ref:sigroteins)", fig.alt="", eval = showSignificant----
sigheat

## ----vennDiagramSig, fig.cap="(ref:vennDiagramSig)", fig.alt="", eval = showSignificant----
xx <- split(ctdata[["ID"]], ctdata$contrast)
if (length(xx) > 1) {
  UpSetR::upset(UpSetR::fromList(xx))
}


## ----sessionInfo--------------------------------------------------------------
pander::pander(sessionInfo())

