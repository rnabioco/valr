#' @details input tbls can be grouped using \code{dplyr::group_by} prior to 
#'   analysis. It is not necessary to group by \code{chrom}because functions do
#'   this internally. Grouping input intervals by \code{strand} will constrain
#'   analyses to the same strand. To compare opposing strands between two tbls,
#'   the strand on the \code{y} tbl can first be be inverted using
#'   \code{flip_strands}.
#'   
#' @seealso \link{flip_strands}
