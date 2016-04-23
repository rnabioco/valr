#' instersect BED intervals
#' 
#' @param x BED intervals 
#' @param y BED intervals 
#' @param max_dist maximum distance between intserections
#' @param strand intersect intervals on same strand
#' 
#' @examples 
#' x <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    500,
#' "chr2", 200,    400,
#' "chr2", 300,    500,
#' "chr2", 800,    900
#' )
#' 
#' y <- tibble::frame_data(
#' ~chrom, ~start, ~end,
#' "chr1", 150,    400,
#' "chr2", 230,    430,
#' "chr2", 350,    430
#' )
#'
#' bed_intersect(x, y)
#' bed_intersect(x, y, max_dist = 50)
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'  
#' @export
bed_intersect <- function(x, y, max_dist = 0, strand = FALSE, strand_opp = FALSE) {
 
  max_rows <- 100
  if (nrow(x) > max_rows || nrow(y) > max_rows) {
    stop("max rows exceeded for r impl of bed_intersect")
  } 
  
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  if (is.null(groups(x)) || groups(x) != "chrom")
    x <- group_by(x, chrom)
  if (is.null(groups(y)) || groups(y) != "chrom")
    y <- group_by(y, chrom)
 
  res <- intersect_rimpl(x, y, max_dist)
  
  if (strand) {
    res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
    res <- filter(res, strand.x != strand.y) 
  }
  
  res
}

#' pure R impl of intersect
intersect_rimpl <- function(x, y, max_dist) {
  full_join(x, y, by = 'chrom') %>%
    rowwise() %>%
    mutate(.overlap = min(end.x, end.y) - max(start.x, start.y)) %>%
    filter(.overlap > max_dist) %>%
    select(-.overlap) %>% ungroup()
}
