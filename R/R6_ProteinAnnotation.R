#' add REV and zz entries - used for testing
#' @param stringsAll character vector of protein IDs
#' @param pattern_decoys prefix for decoy entries
#' @param pattern_contaminants prefix for contaminant entries
#' @export
#'
add_RevCon <- function(
  stringsAll,
  pattern_decoys = "REV_",
  pattern_contaminants = "zz"
) {
  # Set seed for reproducibility
  set.seed(123)

  dd <- data.frame(idx = unique(stringsAll), tomod = unique(stringsAll))
  # Determine the number of elements to prefix
  n <- length(dd$idx)
  num_rev <- ceiling(0.10 * n) # 10 percent with "REV_"
  num_zz <- ceiling(0.05 * n) # 5 percent with "ZZ"

  # Randomly select indices for "REV_" prefix
  indices_rev <- sample(1:n, num_rev, replace = FALSE)

  # Apply the "REV_" prefix
  dd$tomod[indices_rev] <- paste(
    pattern_decoys,
    dd$tomod[indices_rev],
    sep = ""
  )

  # Exclude already modified strings and select indices for "ZZ" prefix
  available_indices <- setdiff(1:n, indices_rev)
  indices_zz <- sample(available_indices, num_zz, replace = FALSE)

  # Apply the "ZZ" prefix
  dd$tomod[indices_zz] <- paste(
    pattern_contaminants,
    dd$tomod[indices_zz],
    sep = ""
  )

  res <- merge(data.frame(idx = stringsAll), dd, by = "idx")
  res <- res[match(stringsAll, res$idx), ] # preserve ordering.
  return(res$tomod)
}

#' simulate peptdata and fitting protein annotation for testing
#' @param Nprot number of proteins to simulate
#' @param PROTEIN if TRUE simulate protein-level data
#' @export
#' @examples
#'
#' res <- sim_data_protAnnot()
#' res <- sim_data_protAnnot(PROTEIN = TRUE)
#'
sim_data_protAnnot <- function(Nprot = 100, PROTEIN = FALSE) {
  if (PROTEIN) {
    istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot)
  } else {
    istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot)
  }
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  tmp_data <- lfqdata$data_long()
  tmp_data$protein_Id <- add_RevCon(tmp_data$protein_Id)
  lfqdata$set_data(tmp_data)
  pids <- grep(
    "^zz|^REV",
    unique(lfqdata$data_long()$protein_Id),
    value = TRUE,
    invert = TRUE
  )
  addannot <- data.frame(
    protein_Id = pids,
    description = stringi::stri_rand_strings(length(pids), 13)
  )
  addannot <- addannot |>
    tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
  # Sample protein lengths from a log-normal distribution
  protein_lengths <- rlnorm(Nprot, meanlog = log(400), sdlog = 0.8)
  protein_lengths <- round(pmax(50, protein_lengths)) # Ensure minimum length of 50 AA
  nr_pep <- round(protein_lengths / 20) # nolint object_usage_linter. example code

  pannot <- ProteinAnnotation$new(
    lfqdata,
    addannot,
    description = "description",
    cleaned_ids = "cleanID",
    pattern_contaminants = "^zz",
    pattern_decoys = "^REV"
  )
  pannot$row_annot$nr_tryptic_peptides <- pannot$row_annot$nr_peptides * 2
  pannot$row_annot$protein_length <- pannot$row_annot$nr_peptides * 10

  return(list(pannot = pannot, lfqdata = lfqdata))
}

#' make lfqdata with row annotation
#' @param Nprot number of proteins to simulate
#' @export
make_annotated_experiment <- function(Nprot = 100) {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  tmp_data <- lfqdata$data_long()
  tmp_data$protein_Id <- add_RevCon(tmp_data$protein_Id)
  lfqdata$set_data(tmp_data)
  pids <- grep(
    "^zz|^REV",
    unique(lfqdata$data_long()$protein_Id),
    value = TRUE,
    invert = TRUE
  )
  addannot <- data.frame(
    protein_Id = pids,
    description = stringi::stri_rand_strings(length(pids), 13)
  )
  addannot <- addannot |>
    tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
  pannot <- ProteinAnnotation$new(
    lfqdata,
    addannot,
    description = "description",
    cleaned_ids = "cleanID",
    pattern_contaminants = "^zz",
    pattern_decoys = "^REV"
  )
  return(list(lfqdata = lfqdata, pannot = pannot))
}


# Decoy / duplicate-ID resolution helpers ----

#' Detect decoy/reverse identifiers (within-duplicate resolution)
#'
#' Thin wrapper delegating to \code{prolfqua::is_decoy} so annotation
#' de-duplication and the quant layer share ONE detector: built-in anchored
#' default prefixes unioned with an optional configured \code{pattern}; an empty
#' / \code{NULL} / no-op (\code{"a^"}) pattern falls back to the defaults only.
#' @param ids character vector of (prefixed) identifiers
#' @param pattern optional configured decoy regex
#' @return logical vector
#' @keywords internal
.detect_decoy_ids <- function(ids, pattern = NULL) {
  prolfqua::is_decoy(ids, pattern = pattern)
}

#' Resolve duplicate protein IDs to one row each
#'
#' Within each duplicated-id group: drop decoy rows when a forward exists (keep
#' the forward), then prefer reviewed \code{sp|} over \code{tr|}, else keep the
#' first. Guarantees one row per id and logs the counts. Standalone decoys (no
#' forward twin) are left untouched.
#' @param row_annot annotation data frame
#' @param pID name of the protein-id column
#' @param full_id name of the column carrying the raw, prefixed id
#' @param pattern_decoys optional configured decoy regex
#' @return \code{row_annot} with one row per \code{pID}
#' @keywords internal
.resolve_unique_protein_ids <- function(row_annot, pID, full_id, pattern_decoys = NULL) {
  ids <- as.character(row_annot[[pID]])
  if (anyDuplicated(ids) == 0L) {
    return(row_annot)
  }
  full <- as.character(row_annot[[full_id]])
  is_decoy <- .detect_decoy_ids(full, pattern_decoys)
  is_sp <- grepl("^sp\\|", full)
  keep <- rep(TRUE, nrow(row_annot))
  dup_ids <- unique(ids[duplicated(ids)])
  n_decoy <- 0L
  n_sp <- 0L
  n_first <- 0L
  for (id in dup_ids) {
    idx <- which(ids == id)
    if (any(!is_decoy[idx])) {
      drop_decoy <- idx[is_decoy[idx]]
      keep[drop_decoy] <- FALSE
      n_decoy <- n_decoy + length(drop_decoy)
      idx <- idx[!is_decoy[idx]]
    }
    if (length(idx) > 1L) {
      sp_idx <- idx[is_sp[idx]]
      if (length(sp_idx) > 0L && length(sp_idx) < length(idx)) {
        keepidx <- sp_idx[1]
        n_sp <- n_sp + 1L
      } else {
        keepidx <- idx[1]
        n_first <- n_first + 1L
      }
      keep[setdiff(idx, keepidx)] <- FALSE
    }
  }
  logger::log_warn(
    "ProteinAnnotation: ",
    length(dup_ids),
    " duplicated '",
    pID,
    "' id(s) collapsed; dropped ",
    n_decoy,
    " decoy row(s); ",
    n_sp,
    " resolved by sp| preference; ",
    n_first,
    " by keep-first."
  )
  row_annot[keep, , drop = FALSE]
}


# ProteinAnnotation ----
#' Decorates LFQData with a row annotation and some protein specific functions.
#'
#' @export
#' @family LFQData
#' @examples
#'
#' istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 100)
#' lfq0 <- prolfqua::LFQData$new(istar$data, istar$config)
#' xd1 <- prolfqua::nr_children_experiment(lfq0$data_long(), lfq0$response(),
#'   lfq0$relevant_hierarchy_keys(), lfq0$file_name(), lfq0$nr_children_col())
#'
#' xd2 <- prolfqua::nr_features_experiment(lfq0$data_long(), lfq0$hierarchy_keys(),
#'   lfq0$relevant_hierarchy_keys())
#' xd1$nr_child_exp |> table()
#'
#' lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
#' tmp <- lfqdata$data_long()
#' tmp$protein_Id <- add_RevCon(tmp$protein_Id)
#' lfqdata$set_data(tmp)
#' pids <- grep("^zz|^REV", unique(lfqdata$data_long()$protein_Id), value = TRUE, invert = TRUE)
#' addannot <- data.frame(
#'   protein_Id = pids,
#'   description = stringi::stri_rand_strings(length(pids), 13)
#' )
#'
#' addannot <- addannot |> tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
#' # ProteinAnnotation$debug("initialize")
#' # debug(nr_obs_sample)
#' xd4 <- prolfqua::nr_obs_sample(lfqdata$data_long(), lfqdata$response(),
#'   lfqdata$relevant_hierarchy_keys(), lfqdata$file_name(), lfqdata$nr_children_col())
#' xd3 <- prolfqua::nr_features_experiment(lfqdata$data_long(), lfqdata$hierarchy_keys(),
#'   lfqdata$relevant_hierarchy_keys())
#'
#' pannot <- ProteinAnnotation$new(lfqdata,
#'   addannot,
#'   description = "description",
#'   cleaned_ids = "cleanID",
#'   pattern_contaminants = "^zz",
#'   pattern_decoys = "^REV"
#' )
#' stopifnot(pannot$annotate_contaminants() == 5)
#' dd <- pannot$clean()
#' pannot$nr_clean()
#' pannot$get_summary()
#' stopifnot(nrow(dd) == 85)
#' tmp <- lfqdata$get_subset(dd)
#' dx2 <- pannot$filter_by_nr_children(exp_nr_children = 2)
#' dx3 <- pannot$filter_by_nr_children(exp_nr_children = 3)
#' stopifnot(nrow(dx2) >= nrow(dx3))
#'
ProteinAnnotation <-
  R6::R6Class(
    "ProteinAnnotation",
    public = list(
      #' @field row_annot data.frame containing further information
      row_annot = NULL,
      #' @field pID column with protein ids
      pID = character(),
      #' @field full_id column with protein id e.g. sp| can be same as pID
      full_id = character(),
      #' @field description name of column containing descriptions
      description = "description",
      #' @field cleaned_ids vector with columns containing addition IDs
      cleaned_ids = character(),
      #' @field exp_nr_children name of columns with the number of peptides
      exp_nr_children = character(),
      #' @field pattern_contaminants pattern_contaminants
      pattern_contaminants = character(),
      #' @field pattern_decoys pattern_decoys
      pattern_decoys = character(),
      #' @description initialize
      #' @param lfqdata data frame from \code{\link[prolfqua]{setup_analysis}}
      #' @param row_annot data frame with row annotation.
      #'   Must have columns matching \code{config$hierarchy_keys_depth()}
      #' @param description name of column with description
      #' @param cleaned_ids names of columns with cleaned Ids
      #' @param full_id column with full protein ID
      #' @param exp_nr_children column with the number of children
      #' @param pattern_contaminants pattern_contaminants
      #' @param pattern_decoys pattern_decoys
      initialize = function(
        lfqdata,
        row_annot = NULL,
        description = NULL,
        cleaned_ids = NULL,
        full_id = NULL,
        exp_nr_children = "nr_peptides",
        pattern_contaminants = NULL,
        pattern_decoys = NULL
      ) {
        self$pID <- lfqdata$relevant_hierarchy_keys()[[1]]
        self$exp_nr_children <- exp_nr_children
        self$pattern_contaminants <- if (is.null(pattern_contaminants)) {
          "a^"
        } else {
          pattern_contaminants
        }
        self$pattern_decoys <- if (is.null(pattern_decoys)) {
          "a^"
        } else {
          pattern_decoys
        }
        self$full_id <- if (!is.null(full_id)) {
          full_id
        } else {
          self$pID
        }
        self$cleaned_ids <- if (!is.null(cleaned_ids)) {
          cleaned_ids
        } else {
          self$pID
        }
        self$description <- if (!is.null(description)) {
          description
        } else {
          self$pID
        }

        self$row_annot <- dplyr::distinct(dplyr::select(lfqdata$data_long(), self$pID))
        if (!is.null(row_annot)) {
          stopifnot(self$pID %in% colnames(row_annot))
          self$row_annot <- dplyr::left_join(
            self$row_annot,
            row_annot,
            by = self$pID
          )
        }
        stopifnot(self$cleaned_ids %in% colnames(self$row_annot))
        stopifnot(self$description %in% colnames(self$row_annot))
        if (!self$exp_nr_children %in% colnames(row_annot)) {
          warning(
            "no exp_nr_children column specified, computing using nr_children_experiment"
          )
          self$row_annot <- dplyr::inner_join(
            self$row_annot,
            prolfqua::nr_children_experiment(
              lfqdata$data_long(),
              response = lfqdata$response(),
              hierarchy_keys_depth = lfqdata$hierarchy_keys()[1],
              file_name = lfqdata$file_name(),
              nr_children_col = lfqdata$nr_children_col(),
              name_nr_child = self$exp_nr_children
            ),
            by = self$pID
          )
        }
        # Invariant: one row per protein ID. Resolve duplicates decoy-aware
        # (drop decoys colliding with a forward; sp| tiebreak; else keep-first).
        self$row_annot <- .resolve_unique_protein_ids(
          self$row_annot,
          self$pID,
          self$full_id,
          self$pattern_decoys
        )
      },
      #' @description
      #' configured decoy pattern, or NULL when none was set
      get_rev_pattern = function() {
        if (
          length(self$pattern_decoys) != 1 ||
            is.na(self$pattern_decoys) ||
            !nzchar(self$pattern_decoys) ||
            identical(self$pattern_decoys, "a^")
        ) {
          return(NULL)
        }
        self$pattern_decoys
      },
      #' @description
      #' annotate contaminants
      #'
      #' Sets the logical \code{CON} column via the shared
      #' \code{prolfqua::is_contaminant} detector (configured pattern unioned with
      #' the built-in defaults; an empty / \code{NULL} / \code{"a^"} pattern falls
      #' back to the defaults only -- never \code{grepl("", x)}, which would flag
      #' every protein). This is the same detector the quant layer uses.
      annotate_contaminants = function() {
        self$row_annot$CON <- prolfqua::is_contaminant(
          as.character(self$row_annot[[self$full_id]]),
          self$pattern_contaminants
        )
        return(sum(self$row_annot$CON))
      },
      #' @description
      #' get summary (contaminants only; decoys are removed at construction)
      get_summary = function() {
        allProt <- nrow(self$row_annot)
        data.frame(
          totalNrOfProteins = allProt,
          percentOfContaminants = round(
            self$annotate_contaminants() / allProt * 100,
            digits = 2
          )
        )
      },
      #' @description number of proteins kept after \code{clean()}
      #' @param contaminants remove contaminants
      nr_clean = function(contaminants = TRUE) {
        nrow(self$clean(contaminants = contaminants))
      },
      #' @description
      #' remove contaminants (always) and, when a decoy pattern was configured,
      #' decoy proteins from the annotation
      #' @param contaminants remove contaminants
      clean = function(contaminants = TRUE) {
        if (contaminants && !("CON" %in% colnames(self$row_annot))) {
          stop("annotate CON")
        }
        res <- self$row_annot
        if (contaminants) {
          res <- res[!res$CON, , drop = FALSE]
        }
        revpat <- self$get_rev_pattern()
        if (!is.null(revpat)) {
          res <- res[
            !grepl(revpat, as.character(res[[self$full_id]])),
            ,
            drop = FALSE
          ]
        }
        res
      },
      #' @description
      #' filter by number children
      #' @param exp_nr_children minimum number of children required
      filter_by_nr_children = function(exp_nr_children = 2) {
        res <- self$row_annot |>
          dplyr::filter(!!sym(self$exp_nr_children) >= exp_nr_children)
        res <- res |> dplyr::select(self$pID, self$exp_nr_children)
        return(res)
      }
    )
  )


#' build Dataset protein annot, defaults are compatible with DIANN
#'
#' @export
#' @param lfqdata LFQData
#' @param msdata data frame
#' @param idcol named vector mapping protein ID column
#' @param cleaned_protein_id column with cleaned protein ID
#' @param protein_description column with protein description
#' @param exp_nr_children column with number of peptides
#' @param full_id column with full protein ID
#' @param more_columns additional columns to include
#' @param pattern_contaminants regex pattern for contaminants
#' @param pattern_decoys regex pattern for decoys
#' @examples
#' # example code
#'
build_protein_annot <- function(
  lfqdata,
  msdata,
  idcol = c("protein_Id" = "Protein.Group"),
  cleaned_protein_id = "Protein.Group.2",
  protein_description = "fasta.header",
  exp_nr_children = "nrPeptides",
  full_id = "fasta.id",
  more_columns = c("fasta.id"),
  pattern_contaminants = "^zz|^CON",
  pattern_decoys = "REV_"
) {
  proteinID_column <- names(idcol)[1]
  msdata <- dplyr::mutate(msdata, !!proteinID_column := !!rlang::sym(idcol))
  length_protIDs <- length(unique(msdata[[proteinID_column]]))
  prot_annot <- dplyr::select(
    msdata,
    dplyr::all_of(unique(c(
      proteinID_column,
      protein_description,
      cleaned_protein_id,
      exp_nr_children,
      full_id,
      more_columns
    )))
  ) |>
    dplyr::distinct()
  stopifnot(length_protIDs == nrow(prot_annot))
  prot_annot <- dplyr::rename(
    prot_annot,
    description = !!rlang::sym(protein_description)
  )
  prot_annot <- dplyr::rename(
    prot_annot,
    IDcolumn = !!rlang::sym(cleaned_protein_id)
  )
  protAnnot <- prolfquapp::ProteinAnnotation$new(
    lfqdata,
    prot_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = full_id,
    exp_nr_children = exp_nr_children,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
  return(protAnnot)
}


#' Dataset protein annot
#'
#' Extracts protein annotation from a data frame, renaming columns and
#' auto-detecting UniProt identifiers. For new code prefer
#' \code{\link{build_protein_annot}} which returns a
#' \code{\link{ProteinAnnotation}} R6 object.
#'
#' @export
#' @param msdata data frame
#' @param idcol named vector mapping protein ID column
#' @param protein_annot fasta header column name
#' @param more_columns more columns to include
dataset_protein_annot <- function(
  msdata,
  idcol = c("protein_Id" = "Protein.Group"),
  protein_annot = "fasta.header",
  more_columns = c("nrPeptides", "fasta.id")
) {
  proteinID_column <- names(idcol)[1]
  msdata <- dplyr::rename(msdata, !!proteinID_column := !!rlang::sym(idcol))
  prot_annot <- dplyr::select(
    msdata,
    dplyr::all_of(c(proteinID_column, protein_annot, more_columns))
  ) |>
    dplyr::distinct()
  prot_annot <- dplyr::rename(
    prot_annot,
    description = !!rlang::sym(protein_annot)
  )

  UNIPROT <- mean(grepl("^sp\\||^tr\\|", prot_annot[[proteinID_column]])) > 0.8
  message("uniprot database : ", UNIPROT)

  if (UNIPROT) {
    prot_annot <- prolfqua::get_uniprot_id_from_fasta_header(
      prot_annot,
      idcolumn = proteinID_column
    )
    prot_annot <- prot_annot |>
      dplyr::rename(!!"IDcolumn" := !!rlang::sym("UniprotID"))
  } else {
    prot_annot$IDcolumn <- prot_annot[[proteinID_column]]
  }
  return(prot_annot)
}
