#' Compute coverage of y intervals over x intervals
#' 
#' @param x tbl of intervals 
#' @param y tbl of intervals 
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strands
#' @param ... extra arguments (not used)
#'
#' @note Book-ended intervals are counted as overlapping by default
#' @family multi-set-ops
#'  
#' @examples 
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' "chr1", 100,    500,
#' "chr2", 200,    400,
#' "chr2", 300,    500,
#' "chr2", 800,    900
#' )
#' 
#' y <- tibble::tribble(
#' ~chrom, ~start, ~end, ~value,
#' "chr1", 150,    400,  100,
#' "chr1", 500,    550,  100,
#' "chr2", 230,    430,  200,
#' "chr2", 350,    430,  300
#' )
#'
#' bed_coverage(x, y)
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html}
#'  
#' @export
bed_coverage <- function(x, y, 
                          strand = FALSE, strand_opp = FALSE, ...) {
  
  if ( ! is_sorted(x) )
    x <- bed_sort(x)
  if ( ! is_sorted(y) )
    y <- bed_sort(y)
  
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)
  
  if (!strand) {
    res <- coverage_impl(x, y)
  } 
  
  if (strand) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- coverage_impl(x_pos, y_pos)
    res_neg <- coverage_impl(x_neg, y_neg)
    res <- bind_rows(res_pos, res_neg)
  } else if (strand_opp) {
    x_pos <- filter(x, strand == "+") 
    y_pos <- filter(y, strand == "+") 
    x_neg <- filter(x, strand == "-") 
    y_neg <- filter(y, strand == "-") 
    res_pos <- coverage_impl(x_pos, y_neg)
    res_neg <- coverage_impl(x_neg, y_pos)
    res <- bind_rows(res_pos, res_neg)
  }
  
  res
}
