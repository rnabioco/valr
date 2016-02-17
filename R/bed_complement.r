#' identify intervals in a genome that are not covered by a query
#' 
#' @param bed_df BED data.frame
#' @param genome chrom sizes
#' 
#' @return \code{data.frame}
#' 
#' @examples 
#' 
#' bed_df <- dplyr::tibble(
#'    ~chrom, ~start, ~end,
#'    "chr1", 100,    300,
#'    "chr2", 200,    400,
#'    "chr3", 500,    600
#' )
#' 
#' bed_complement(bed_df)
#' 
#' @export
bed_complement <- function(bed_df, genome) {
  res <- bed_df %>% 
    group_by(chrom) %>%
    bed_merge() %>%
    bed_complement_(., genome) %>%
    ungroup()
}

bed_complement_ <- function(df, genome) {
  res <- complement_cpp(df, genome)
  res
}
