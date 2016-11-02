#' @details input tbls can be grouped using \code{dplyr::group_by} prior to 
#'   analysis. Grouping by \code{chrom} is done by default. Grouping by 
#'   \code{strand} will constrain analyses to the same strand. To compare 
#'   opposing strands between two tbls, strands on the \code{y} tbl can first be
#'   be inverted using \code{flip_strands}.
#'   
#' @seealso \link{flip_strands}
