#' read annotation files
#' @return list with annot (annotation table), atable (analtysis table configuration), contrasts list with contrasts.
#' @param dsf annotation table
#' @param repeated is this a repeated measurement
#' @param SAINT is this a SAINTexpress analysis
#' @export
read_annotation <- function(dsf, repeated = TRUE, SAINT = FALSE, prefix = "G_", QC = FALSE){
 res <- AnnotationProcessor$new(repeated = repeated, SAINT = SAINT, prefix = prefix,QC = QC)$read_annotation(dsf)
 return(res)
}

#' extract contrast from annotation file
#' @export
#' @examples
#'
#' annot <- data.frame(names = c("a1","b1"), group= c("a","b"), ddd = c("T","C"))
#' testthat::expect_error(extract_contrasts(annot))
#' annot$control <- annot$ddd
#' contrast <- extract_contrasts(annot)
#' stopifnot(contrast == "G_a - G_b")
#'
#' annot$Contrast <- c("G_a - G_b","G_b - G_a")
#' annot$ContrastName <- c("a_vs_b","b_vs_a")
#' annot$control <- NULL
#' ct <- extract_contrasts(annot)
#' stopifnot(length(ct) == 2)
extract_contrasts <- function(annot, prefix = "G_", group = "group") {
  AnnotationProcessor$new(prefix = prefix)$extract_contrasts(annot, group = group)
}

#' add vector of contrasts to annotation data frame
#' @export
#' @examples
#' annot <- data.frame(Group = rep(c("A","B","C"), each = 3))
#' annot$Name
add_contrasts_vec <- function(xx, Contrasts){
  AnnotationProcessor$new()$add_contrasts_vec(xx, Contrasts)
}

