#' identify intervals in a genome that are not covered by a query
#' 
#' @param bed_tbl tbl of intervals
#' @param genome chrom sizes
#' 
#' @return \code{data.frame}
#' 
#' @examples 
#' genome <- dplyr::tibble(
#'    ~chrom,  ~size,
#'    "chr1", 500,
#'    "chr2", 600,
#'    "chr3", 800
#' ) 
#' 
#' bed_tbl <- dplyr::tibble(
#'    ~chrom, ~start, ~end,
#'    "chr1", 100,    300,
#'    "chr1", 200,    400,
#'    "chr2", 1,      100,
#'    "chr2", 200,    400,
#'    "chr3", 500,    600
#' )
#' 
#' # intervals not covered by bed_tbl
#' bed_complement(bed_tbl, genome)
#' 
#' @export
bed_complement <- function(bed_tbl, genome) {

  if ( ! is_merged(bed_tbl) ) {
    res <- bed_merge(bed_tbl) %>% mutate(chrom = as.character(chrom))
  } 

  res <- complement_impl(res, genome) 
  
  res <- bed_sort(res)
 
  res 
}
