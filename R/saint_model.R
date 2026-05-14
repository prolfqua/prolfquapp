.saint_first_existing <- function(data, candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates[candidates %in% colnames(data)][1]
}

.saint_bait_col <- function(lfq_data) {
  data <- lfq_data$data_long()
  factor_cols <- lfq_data$relevant_factor_keys()
  bait_col <- .saint_first_existing(
    data,
    c(
      factor_cols[grepl("^bait", factor_cols, ignore.case = TRUE)],
      "bait",
      "Bait"
    )
  )
  if (is.na(bait_col)) {
    stop(
      "SAINT model requires a bait column. Read annotation with SAINT = TRUE ",
      "or provide an annotation column starting with 'bait'.",
      call. = FALSE
    )
  }
  bait_col
}

.saint_control_col <- function(data) {
  control_col <- .saint_first_existing(
    data,
    c(
      "CONTROL",
      grep("^control", colnames(data), value = TRUE, ignore.case = TRUE)
    )
  )
  if (is.na(control_col)) {
    stop(
      "SAINT model requires a CONTROL column with C/T values.",
      call. = FALSE
    )
  }
  invalid <- setdiff(unique(stats::na.omit(data[[control_col]])), c("C", "T"))
  if (length(invalid) > 0) {
    stop(
      "SAINT CONTROL column must contain only C and T values.",
      call. = FALSE
    )
  }
  control_col
}

.saint_prepare_data <- function(lfq_data, row_annot) {
  data <- lfq_data$data_long()
  protein_id <- lfq_data$relevant_hierarchy_keys()[[1]]
  response <- lfq_data$response()

  if (!protein_id %in% colnames(data)) {
    stop("SAINT model requires a protein hierarchy column.", call. = FALSE)
  }
  if (!response %in% colnames(data)) {
    stop("SAINT model response column not found: ", response, call. = FALSE)
  }

  annot_cols <- setdiff(colnames(row_annot), colnames(data))
  if (length(annot_cols) > 0) {
    annot <- dplyr::distinct(row_annot[,
      c(protein_id, annot_cols),
      drop = FALSE
    ])
    data <- dplyr::left_join(data, annot, by = protein_id, multiple = "all")
  }

  length_col <- .saint_first_existing(
    data,
    c("protein.length", "protein_length", "Protein.Length", "proteinLength")
  )
  if (is.na(length_col)) {
    data$.saint_protein_length <- 400L
  } else {
    data$.saint_protein_length <- suppressWarnings(as.integer(data[[
      length_col
    ]]))
    fallback <- mean(data$.saint_protein_length, na.rm = TRUE)
    if (is.nan(fallback)) {
      fallback <- 400L
    }
    data$.saint_protein_length[is.na(data$.saint_protein_length)] <- as.integer(
      fallback
    )
  }

  gene_col <- .saint_first_existing(
    data,
    c("geneNames", "Gene.names", "Gene", "gene", "GN", "cleanID", "description")
  )
  if (is.na(gene_col)) {
    data$.saint_gene_names <- data[[protein_id]]
  } else {
    data$.saint_gene_names <- data[[gene_col]]
    data$.saint_gene_names[is.na(data$.saint_gene_names)] <- data[[
      protein_id
    ]][is.na(data$.saint_gene_names)]
  }

  list(
    data = data,
    protein_id = protein_id,
    response = response,
    sample_col = lfq_data$file_name(),
    bait_col = .saint_bait_col(lfq_data),
    control_col = .saint_control_col(data)
  )
}

.build_saint_contrast_result <- function(
  lfq_data,
  row_annot,
  spc = FALSE,
  engine = "r"
) {
  prepared <- .saint_prepare_data(lfq_data, row_annot$row_annot)
  saint_input <- prolfquasaint::protein_2localSaint(
    prepared$data,
    quantcolumn = prepared$response,
    proteinID = prepared$protein_id,
    geneNames = ".saint_gene_names",
    proteinLength = ".saint_protein_length",
    IP_name = prepared$sample_col,
    baitCol = prepared$bait_col,
    CorTCol = prepared$control_col
  )
  saint_result <- prolfquasaint::runSaint(
    saint_input,
    spc = spc,
    engine = engine,
    CLEANUP = TRUE
  )
  saint_result$list[[prepared$protein_id]] <- saint_result$list$Prey
  contrast <- prolfquasaint::ContrastsSAINTexpress$new(
    saint_result$list,
    subject_id = prepared$protein_id
  )
  list(
    contrast = contrast,
    input = saint_input,
    result = saint_result,
    protein_id = prepared$protein_id
  )
}
