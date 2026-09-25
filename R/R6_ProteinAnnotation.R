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
  set.seed(123)
  ids <- unique(stringsAll)
  n <- length(ids)
  # 10 percent get the decoy prefix, another 5 percent the contaminant prefix
  indices_rev <- sample(1:n, ceiling(0.10 * n), replace = FALSE)
  indices_zz <- sample(setdiff(1:n, indices_rev), ceiling(0.05 * n), replace = FALSE)
  tomod <- ids
  tomod[indices_rev] <- paste0(pattern_decoys, tomod[indices_rev])
  tomod[indices_zz] <- paste0(pattern_contaminants, tomod[indices_zz])
  tomod[match(stringsAll, ids)]
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
  istar <- if (PROTEIN) {
    prolfqua::sim_lfq_data_protein_config(Nprot = Nprot)
  } else {
    prolfqua::sim_lfq_data_peptide_config(Nprot = Nprot)
  }
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  tmp_data <- lfqdata$data_long()
  tmp_data$protein_Id <- add_RevCon(tmp_data$protein_Id)
  lfqdata$set_data(tmp_data)
  pids <- grep("^zz|^REV", unique(lfqdata$data_long()$protein_Id), value = TRUE, invert = TRUE)
  addannot <- data.frame(protein_Id = pids, description = stringi::stri_rand_strings(length(pids), 13)) |>
    tidyr::separate(protein_Id, c("cleanID", NA), remove = FALSE)
  pannot <- ProteinAnnotation$new(
    lfqdata,
    addannot,
    description = "description",
    cleaned_ids = "cleanID",
    pattern_contaminants = "^zz",
    pattern_decoys = "^REV"
  )
  pannot$row_annot$nr_tryptic_peptides <- pannot$row_annot$nrPeptides * 2
  pannot$row_annot$protein_length <- pannot$row_annot$nrPeptides * 10
  list(pannot = pannot, lfqdata = lfqdata)
}

# Decoy / duplicate-ID resolution helpers ----

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
#' @noRd
.resolve_unique_protein_ids <- function(row_annot, pID, full_id, pattern_decoys = NULL) {
  # pID may name several columns -- protein_Id and site for a site-level
  # analysis -- in which case uniqueness is a property of the combination, not
  # of the protein: one protein legitimately carries many sites.
  ids <- if (length(pID) == 1L) {
    as.character(row_annot[[pID]])
  } else {
    do.call(paste, c(lapply(pID, function(k) as.character(row_annot[[k]])), sep = "\r"))
  }
  if (anyDuplicated(ids) == 0L) {
    return(row_annot)
  }
  full <- as.character(row_annot[[full_id]])
  is_decoy <- prolfqua::is_decoy(full, pattern = pattern_decoys)
  is_sp <- grepl("^sp\\|", full)
  keep <- rep(TRUE, nrow(row_annot))
  dup_ids <- unique(ids[duplicated(ids)])
  n_decoy <- 0L
  n_sp <- 0L
  n_first <- 0L
  for (id in dup_ids) {
    idx <- which(ids == id)
    if (any(!is_decoy[idx])) {
      keep[idx[is_decoy[idx]]] <- FALSE
      n_decoy <- n_decoy + sum(is_decoy[idx])
      idx <- idx[!is_decoy[idx]]
    }
    if (length(idx) > 1L) {
      sp_idx <- idx[is_sp[idx]]
      by_sp <- length(sp_idx) > 0L && length(sp_idx) < length(idx)
      n_sp <- n_sp + by_sp
      n_first <- n_first + !by_sp
      keep[setdiff(idx, if (by_sp) sp_idx[1] else idx[1])] <- FALSE
    }
  }
  logger::log_warn(
    "ProteinAnnotation: {length(dup_ids)} duplicated '{paste(pID, collapse = ' + ')}' id(s) collapsed; ",
    "dropped {n_decoy} decoy row(s); {n_sp} resolved by sp| preference; {n_first} by keep-first."
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
      #' @field pID key column(s) of the annotation: the protein id, plus the
      #'   site id when the annotation describes sites
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
      #' @param row_annot data frame with row annotation. Must carry the
      #'   protein-id column; when it also carries a deeper hierarchy key such
      #'   as \code{site}, the annotation is taken to describe rows at that
      #'   level and stays one row per key combination.
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
        exp_nr_children = "nrPeptides",
        pattern_contaminants = NULL,
        pattern_decoys = NULL
      ) {
        # The annotation is a row annotation keyed on whatever identifies a row
        # of the analysis: protein_Id, or protein_Id and site together when the
        # supplied table also carries the site key. Everything that is a property
        # of the protein alone keys off the first element.
        keys <- lfqdata$relevant_hierarchy_keys()
        self$pID <- keys[[1]]
        if (!is.null(row_annot)) {
          if (!keys[[1]] %in% colnames(row_annot)) {
            stop(
              "row_annot must carry the '",
              keys[[1]],
              "' column; it has: ",
              paste(colnames(row_annot), collapse = ", ")
            )
          }
          self$pID <- keys[keys %in% colnames(row_annot)]
        }
        self$exp_nr_children <- exp_nr_children
        self$pattern_contaminants <- if (is.null(pattern_contaminants)) "a^" else pattern_contaminants
        self$pattern_decoys <- if (is.null(pattern_decoys)) "a^" else pattern_decoys
        self$full_id <- if (is.null(full_id)) self$pID[[1]] else full_id
        self$cleaned_ids <- if (is.null(cleaned_ids)) self$pID[[1]] else cleaned_ids
        self$description <- if (is.null(description)) self$pID[[1]] else description

        self$row_annot <- dplyr::distinct(dplyr::select(lfqdata$data_long(), dplyr::all_of(self$pID)))
        if (!is.null(row_annot)) {
          self$row_annot <- dplyr::left_join(self$row_annot, row_annot, by = self$pID)
        }
        stopifnot(self$cleaned_ids %in% colnames(self$row_annot))
        stopifnot(self$description %in% colnames(self$row_annot))
        if (!self$exp_nr_children %in% colnames(row_annot)) {
          warning("no exp_nr_children column specified, computing using nr_children_experiment")
          self$row_annot <- dplyr::inner_join(
            self$row_annot,
            prolfqua::nr_children_experiment(
              lfqdata$data_long(),
              response = lfqdata$response(),
              hierarchy_keys_depth = self$pID,
              file_name = lfqdata$file_name(),
              nr_children_col = lfqdata$nr_children_col(),
              name_nr_child = self$exp_nr_children
            ),
            by = self$pID
          )
        }
        # Invariant: one row per protein ID. Resolve duplicates decoy-aware
        # (drop decoys colliding with a forward; sp| tiebreak; else keep-first).
        self$row_annot <- .resolve_unique_protein_ids(self$row_annot, self$pID, self$full_id, self$pattern_decoys)
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
      #' Sets the logical \code{CON} column via \code{prolfqua::is_contaminant},
      #' the same detector the quant layer uses.
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
        res <- self$row_annot
        if (contaminants) {
          if (!"CON" %in% colnames(res)) {
            stop("annotate CON")
          }
          res <- res[!res$CON, , drop = FALSE]
        }
        revpat <- self$get_rev_pattern()
        if (!is.null(revpat)) {
          res <- res[!grepl(revpat, as.character(res[[self$full_id]])), , drop = FALSE]
        }
        res
      },
      #' @description
      #' filter by number children
      #' @param exp_nr_children minimum number of children required
      filter_by_nr_children = function(exp_nr_children = 2) {
        self$row_annot |>
          dplyr::filter(!!sym(self$exp_nr_children) >= exp_nr_children) |>
          dplyr::select(dplyr::all_of(c(self$pID, self$exp_nr_children)))
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
    description = !!rlang::sym(protein_description),
    IDcolumn = !!rlang::sym(cleaned_protein_id)
  )
  ProteinAnnotation$new(
    lfqdata,
    prot_annot,
    description = "description",
    cleaned_ids = "IDcolumn",
    full_id = full_id,
    exp_nr_children = exp_nr_children,
    pattern_contaminants = pattern_contaminants,
    pattern_decoys = pattern_decoys
  )
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
  prot_annot <- dplyr::rename(prot_annot, description = !!rlang::sym(protein_annot))
  UNIPROT <- mean(grepl("^sp\\||^tr\\|", prot_annot[[proteinID_column]])) > 0.8
  message("uniprot database : ", UNIPROT)

  if (UNIPROT) {
    prot_annot <- prolfqua::get_uniprot_id_from_fasta_header(prot_annot, idcolumn = proteinID_column) |>
      dplyr::rename(IDcolumn = "UniprotID")
  } else {
    prot_annot$IDcolumn <- prot_annot[[proteinID_column]]
  }
  prot_annot
}
