#' @details input tbls are grouped by \code{chrom} by default, and additional
#'   groups can be added using \code{\link[dplyr]{group_by}}. For example,
#'   grouping by \code{strand} will constrain analyses to the same strand. To
#'   compare opposing strands across two tbls, strands on the \code{y} tbl can
#'   first be inverted using \code{\link{flip_strands}}.
