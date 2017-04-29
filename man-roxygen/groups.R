#' @md
#' @details input tbls are grouped by `chrom` by default, and additional
#'   groups can be added using [dplyr::group_by()]. For example,
#'   grouping by `strand` will constrain analyses to the same strand. To
#'   compare opposing strands across two tbls, strands on the `y` tbl can
#'   first be inverted using [flip_strands()].
