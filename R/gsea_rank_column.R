# The column to rank a GSEA file on: the backend's test statistic, which is
# signed and unbounded. A backend that reports no p-value scores a bounded
# probability instead (SAINTexpress: SaintScore), which carries no direction and
# cannot order a ranked list, so NULL asks it for its own effect size.
.gsea_rank_column <- function(cfg) {
  if (isTRUE(cfg$has_pvalue())) cfg$score_col else NULL
}
