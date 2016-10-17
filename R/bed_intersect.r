#' Identify intersecting intervals.
#' 
#' @param x tbl of intervals 
#' @param y tbl of intervals 
#' @param invert report `x` intervals not in `y`
#' @param strand intersect intervals on same strand
#' @param strand_opp intersect intervals on opposite strands
#' @param suffix colname suffixes in output
#' @param ... extra arguments (not used)
#'
#' @note Book-ended intervals have \code{.overlap} values of 0 in the output.
#'  
#' @examples 
#' 
#' x <- tibble::tribble(
#' ~chrom, ~start, ~end,
#' 'chr1', 25,      50,
#' 'chr1', 100,     125
#' )
#' 
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75
#' )
#' 
#' bed_glyph(bed_intersect(x, y))
#' 
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
#' bed_intersect(x, y)
#'  
#' @family multi-set-ops
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'  
#' @export
bed_intersect <- function(x, y, invert = FALSE,
                          strand = FALSE, strand_opp = FALSE,
                          suffix = c('.x', '.y'), ...) {
  
  # check_suffix
  if (!is.character(suffix) || length(suffix) != 2)
    stop("`suffix` must be a character vector of length 2.", call. = FALSE)
 
  if (strand && !('strand' %in% colnames(x) && 'strand' %in% colnames(y)))
    stop("`strand` specified on unstranded data_frame", call. = FALSE)
  
 # if ( ! is_sorted(x) )
#    x <- bed_sort(x)
 # if ( ! is_sorted(y) )
#    y <- bed_sort(y)
 
  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)

  suffix <- list(x = suffix[1], y = suffix[2])

  res <- intersect_impl(x, y, suffix$x, suffix$y)
  
  # XXX: probably better to group matched / mismatched strands and intersect
  # among those
  if (strand) {
     res <- filter(res, strand.x == strand.y) 
  } else if (strand_opp) {
     res <- filter(res, strand.x != strand.y) 
  }
 
  if (invert) {
    colspec <- c('chrom' = 'chrom', 'start' = 'start.x', 'end' = 'end.x') 
    res <- anti_join(x, res, by = colspec)
    res <- ungroup(res)
  }
  
  res
}
