#' Sort intervals
#' 
#' @param intervals tbl of intervals
#' 
#' @examples 
#' bed_df <- dplyr::tibble(
#'    ~chrom, ~start, ~end,
#'    "chr8", 500, 1000,
#'    "chr8", 1000, 5000,
#'    "chr1", 100, 300,
#'    "chr1", 100, 200
#' )
#' 
#' bed_df %>% bed_sort()
#' @export
bed_sort <- function(intervals) {
   res <- intervals %>%
     arrange(chrom, start, end)
   res
}