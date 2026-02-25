params <-
list(configuration = prolfqua::prolfqua_data("data_ionstar")$filtered()$config$clone(deep = TRUE), 
    data = prolfqua::prolfqua_data("data_ionstar")$filtered()$data, 
    plot_density = TRUE, plot_sd_vs_mean = FALSE, pep = TRUE)

## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  fig.width = 5.5,
  fig.height = 3.5,
  fig.align = "center",
  fig.pos = "H"
)

#params <- prr
data <- eval(params$data)

local_project_conf <- params$project_conf
local_config <- eval(params$configuration)
if (!is.null(local_config$is_intensity_transformed)) {
  local_config <- prolfqua::old2new(local_config)
}




if (is.null(local_project_conf)) {
  local_project_conf <- list()
}

if (is.null(params$plot_density)) {
  params$plot_density <- TRUE
}

if (is.null(params$plot_sd_vs_mean)) {
  params$plot_sd_vs_mean <- FALSE
}
old.theme <- ggplot2::theme_set(ggplot2::theme_classic())


## ----hierarchyCounts----------------------------------------------------------
library(rlang)
library(prolfqua)
lfqdata <- prolfqua::LFQData$new(data, local_config)

x <- lfqdata$hierarchy_counts()

prolfqua::table_facade( 
  data.frame(NR = x),
  caption = paste0("Nr of ",
                   if(params$pep){"proteins and peptides"}
                   else
                   {"proteins"} ,
                   " detected in all samples."))


## ----defineNrSamples----------------------------------------------------------

nrSamples <- nrow(lfqdata$factors())
width = max(8, nrSamples*0.2)

## ----hierarchyCountsSampleBarplot, fig.cap="(ref:hierarchyCountsSampleBarplot)", fig.width=width, fig.height=6----
summarizer <- lfqdata$get_Summariser()
summarizer$plot_hierarchy_counts_sample()


## ----summarizeHier------------------------------------------------------------
x3 <- lfqdata$summarize_hierarchy()


## ----defineHeight-------------------------------------------------------------
height <- max(7, nrSamples * 0.1)

## ----missingFigIntensityHistorgram, fig.width=height, fig.height=height, fig.cap="(ref:missingFigIntensityHistorgram)"----

plotter <- lfqdata$get_Plotter()
p <- plotter$missigness_histogram()
xx3 <- summarizer$plot_missingness_per_group()
xx4 <- summarizer$plot_missingness_per_group_cumsum()
gridExtra::grid.arrange(p, xx3, xx4, ncol = 1)


## ----preparemissingnessHeatmap, include=FALSE---------------------------------
pNAmap <- plotter$NA_heatmap()
width <- max(8, nrSamples * 0.15)

## ----missingnessHeatmap, fig.width=width, fig.height=width, fig.align='center', fig.cap="(ref:missingnessHeatmap)"----
print(pNAmap)

## ----checktransformation, include = FALSE-------------------------------------
show_text <- !lfqdata$is_transformed()

## ----variability, eval = show_text, results='asis'----------------------------
cat("## Variablity of the raw intensities\n\n")

cat("Without applying any intensity scaling and data preprocessing, the peptide intensities in all samples should be similar. To asses this we plotted the distribution of the peptide intensities in the samples (Figure \\@ref(fig:plotDistributions)) as well as the distribution of the coefficient of variation CV for all peptides in the samples (Figure \\@ref(fig:intensityDistribution)). Table \\@ref(tab:printTable) summarises the CV.")


## ----plotDistributions, fig.cap="Density plot of peptide level Coefficient of Variations (CV).", fig.height=6, fig.width=8 , eval = show_text----
stats <- lfqdata$get_Stats()
if (params$plot_density) {
  p1 <- stats$density() +
    ggplot2::labs(tag = "A") +
    ggplot2::xlim(0, 150) #+
    #ggplot2::theme(legend.position = "none")
  p2 <- stats$density_median() +
    ggplot2::labs(tag = "B") +
    ggplot2::xlim(0, 150) +
    ggplot2::theme(legend.position = "bottom")
  gridExtra::grid.arrange(p1,p2)
} else {
  p1 <- stats$violin() + ggplot2::labs(tag='A')
  p1
  p2 <- stats$violin_median() + ggplot2::labs(tag='B')
  gridExtra::grid.arrange(p1,p2)
}

## ----eval=show_text-----------------------------------------------------------
if (params$plot_sd_vs_mean) {
  p0 <- stats$stdv_vs_mean() + labs(tag='A')
}

## ----computeCVQuantiles, include=FALSE, eval = show_text----------------------
cv_quantiles_res <- stats$stats_quantiles(probs = c(0.5, 0.6, 0.7, 0.8, 0.9))$wide

## ----printTable, results="asis", eval = show_text-----------------------------
prolfqua::table_facade(cv_quantiles_res,
                         caption = "Summary of the coefficient of variation (CV) at the 50th, 60th, 70th, 80th and 90th percentile.")



## ----intensityDistribution, fig.cap="Distribution of unnormalized intensities.",fig.height = 6, fig.width=width, eval = show_text----
if(params$plot_density){
  p0 <- plotter$intensity_distribution_density() + 
    ggplot2::theme(legend.text = ggplot2::element_text(size=5))
  p0
} else{
  p1 <- plotter$intensity_distribution_violin()
  p1
}


## ----transformIntensities-----------------------------------------------------
if (show_text) {
  tr <- lfqdata$get_Transformer()
  tr <- tr$log2()
  tr <- tr$robscale()
  dataTransformed <- tr$lfq
} else{
  dataTransformed <- lfqdata
}


## ----plotTransformedIntensityDistributions, fig.cap="(ref:plotTransformedIntensityDistributions)", fig.height = 6, fig.width= width----
plotter <- dataTransformed$get_Plotter()

if (params$plot_density) {
  plotter$intensity_distribution_density() + 
    ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
} else {
  plotter$intensity_distribution_violin()
}

## ----preparecorrelationHeat, include= FALSE-----------------------------------
chmap <- plotter$heatmap_cor()
width <- max(8, nrSamples * 0.15)

## ----correlationHeat, fig.height = 5, fig.cap="(ref:correlationHeat)", fig.width=width, fig.height=width----
print(chmap)

## ----pairsplotSmooth, fig.cap = "Pairsplot - scatterplot of samples.", fig.height=12, fig.width=12----
plotter$pairs_smooth(max = 10)

## ----sdviolinplots,fig.cap="(ref:sdviolinplots)", fig.height=6, fig.width=8----

st <- dataTransformed$get_Stats()

if (params$plot_density) {
  p1 <- st$density()
    ggplot2::labs(tag = "A") +
    ggplot2::theme(legend.position = "none")
  p2 <-
    st$density_median() +
    ggplot2::labs(tag = "B") +
    ggplot2::theme(legend.position = "bottom")
  gridExtra::grid.arrange(p1, p2)
} else {
  p1 <- st$violin() + 
    ggplot2::labs(tag = "A")
  p2 <- st$violin_median()  + 
    ggplot2::labs(tag = "B")
  gridExtra::grid.arrange(p1, p2)
}

## ----sdecdf,fig.cap="(ref:sdecdf)", fig.height=6, fig.width=8-----------------
  p1 <- st$density(ggstat = "ecdf")
    ggplot2::labs(tag = "A") +
    ggplot2::theme(legend.position = "none")
  p2 <-
    st$density_median(ggstat = "ecdf") +
    ggplot2::theme(legend.position = "bottom")
  gridExtra::grid.arrange(p1, p2)

## ----fig.cap="Standard Deviation vs Mean"-------------------------------------
if (params$plot_sd_vs_mean) {
  st$stdv_vs_mean()
}

## ----computeSDQuantiles, include=FALSE----------------------------------------
sd_quantile_res2 <- st$stats_quantiles(probs = c(0.5, 0.6, 0.7, 0.8, 0.9))$wide

## ----printSDTable-------------------------------------------------------------
prolfqua::table_facade(sd_quantile_res2,
  caption = "Summary of the distribution of standard deviations at the 50th, 60th, 70th, 80th and 90th percentile.")


## ----prepareoverviewHeat, include = FALSE-------------------------------------
hm <- plotter$heatmap()
width = max(10, nrSamples * 0.15)

## ----overviewHeat, fig.cap="(ref:overviewHeat)", fig.height=10, fig.width=width----
print(hm)

## ----figSampleSize, fig.cap="(ref:figSampleSize)", fig.width=8, fig.alt=""----

sampleSize2 <- st$power_t_test_quantiles(probs = c(0.5,0.75), delta = c(0.59,1,2))
nudgeval <- max(sampleSize2$N) * 0.05

sampleSize2 <- sampleSize2 |> dplyr::mutate(delta = paste0("delta = ", delta))

#stats_res <- summarize_stats(dataTransformed$data, dataTransformed$config)
#xx <- summarize_stats_quantiles(stats_res, local_config,probs = c(0.5,0.75))
#lfq_power_t_test_quantiles_V2(xx$long,delta = c(0.59,1,2))

ggplot2::ggplot(sampleSize2, ggplot2::aes(x = probs, y = N)) +
  ggplot2::geom_bar(stat = "identity", color = "black", fill = "white") +
  ggplot2::geom_text(ggplot2::aes(label = N), nudge_y = nudgeval) + 
  ggplot2::facet_wrap(as.formula(paste0("~ ", "delta + " , paste0( lfqdata$config$factor_keys_depth(), collapse = " + " ))))
   
tmp <- sampleSize2 |> tidyr::pivot_wider(id_cols = c("probs","sdtrimmed", lfqdata$config$factor_keys_depth()), names_from =  delta, values_from = N  )


## ----sampleSize---------------------------------------------------------------
prolfqua::table_facade(tmp,
  caption = "Sample size needed to detect a difference log fold change greater than delta with a significance level of 0.05 and power 0.8 when using a t-test to compare means.")


## ----sampleNameRawFileMapping-------------------------------------------------
prolfqua::table_facade(
  dataTransformed$factors(),
  caption = "Mapping of raw file names to sample names used throughout this report.")


## ----hierarchyCountsSample----------------------------------------------------
caption <- paste("Number of quantified ",
                 if (params$pep) { "peptides and proteins" } else {"peptides" }, " per sample.")

st <- dataTransformed$get_Summariser()
prolfqua::table_facade(
  st$hierarchy_counts_sample()
  , caption = caption)


## ----resetTheme, include=FALSE------------------------------------------------
ggplot2::theme_set(old.theme)

