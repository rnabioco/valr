#' @details input tbls can be grouped using \code{\link[dplyr]{group_by}} prior to 
#'   analysis. Input tbls are grouped by \code{chrom} by default. Grouping by 
#'   \code{strand} will constrain analyses to the same strand. To compare 
#'   opposing strands between two tbls, strands on the \code{y} tbl can first be
#'   be inverted using \code{\link{flip_strands}}.
