# Shared fixture: a minimal DEA result SummarizedExperiment shaped exactly like
# the one DEAReportGenerator$make_SummarizedExperiment() writes -- an annotation
# rowData frame, results-only contrast frames, and the metadata schema.

make_dea_summarized_experiment <- function() {
  feature_names <- c("P1~S10", "P2~S20")
  sample_names <- c("S1", "S2", "S3")
  raw <- matrix(
    c(10, NA, 30, 20, 40, 60),
    nrow = 2,
    dimnames = list(feature_names, sample_names)
  )
  transformed <- log2(raw)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rawData = raw,
      transformedData = transformed,
      nr_children = matrix(
        c(2, 3, 2, 3, 2, 3),
        nrow = 2,
        dimnames = list(feature_names, sample_names)
      )
    ),
    colData = S4Vectors::DataFrame(
      sampleName = sample_names,
      group = c("A", "A", "B"),
      row.names = sample_names
    ),
    metadata = list(
      artifact_type = "dea_results",
      schema_version = "2.0.0",
      source_software = "DIANN",
      feature_keys = c("protein_Id", "site"),
      sample_key = "sampleName",
      bfabric_urls = list(projectURL = "https://example.org/project/1"),
      provenance = list(software = "DIANN", workunit_Id = 42),
      contrasts = data.frame(
        contrast_name = c("A/B", "A%2FB"),
        contrast = c("A - B", "A + B")
      ),
      formula = data.frame(formula = "abundance ~ group"),
      default_model = "lm",
      analysis_configuration_raw = list(
        sample_name = "sampleName",
        hierarchy = list(protein_Id = "protein", site = "site")
      ),
      analysis_configuration = list(
        sample_name = "sampleName",
        hierarchy = list(protein_Id = "protein", site = "site")
      ),
      contrast_configuration = list(
        subject_id = "protein_Id",
        model_name_col = "modelName",
        contrast_col = "contrast",
        effect_col = "diff",
        score_col = "statistic",
        pvalue_col = "p.value",
        fdr_col = "FDR",
        avg_abundance_col = "avgAbd"
      )
    )
  )
  annotation <- data.frame(
    protein_Id = c("P1", "P2"),
    site = c("S10", "S20"),
    SequenceWindow = c("AAAAASAAAAA", "BBBBBSBBBBB"),
    row.names = feature_names
  )
  contrast_result <- function(contrast, diff) {
    data.frame(
      modelName = "lm",
      estimate_type = "observed",
      contrast = contrast,
      diff = diff,
      statistic = diff * 2,
      p.value = c(0.01, 0.2),
      FDR = c(0.02, 0.2),
      row.names = feature_names
    )
  }
  SummarizedExperiment::rowData(se)[["annotation"]] <- annotation
  SummarizedExperiment::rowData(se)[["constrast_A/B"]] <-
    contrast_result("A/B", c(1, -1))
  SummarizedExperiment::rowData(se)[["constrast_A%2FB"]] <-
    contrast_result("A%2FB", c(2, -2))
  SummarizedExperiment::rowData(se)[["stats_normalized_wide"]] <-
    data.frame(
      protein_Id = c("P1", "P2"),
      site = c("S10", "S20"),
      meanAbundance_All = c(5, 6),
      row.names = feature_names
    )
  SummarizedExperiment::rowData(se)[["stats_raw_wide"]] <-
    data.frame(
      protein_Id = c("P1", "P2"),
      site = c("S10", "S20"),
      meanAbundance_All = c(25, 40),
      row.names = feature_names
    )
  se
}
