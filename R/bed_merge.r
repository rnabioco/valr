#' Merge overlapping intervals.
#' 
#' Operations can be performed on merged intervals by specifying name-value
#' pairs.
#' 
#' @param x tbl of intervals
#' @param max_dist maximum distance between intervals to merge
#' @param ... name-value pairs that specify operations on merged intervals
#'   
#' @template groups
#'   
#' @return \code{data_frame}
#' @family single-set-ops
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#'   
#' @examples 
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',  1,      50,
#'   'chr1',  10,     75,
#'   'chr1',  100,    120
#' )
#' 
#' bed_glyph(bed_merge(x))
#' 
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end, ~value, ~strand,
#'  "chr1", 1,      50,   1,      '+',
#'  "chr1", 100,    200,  2,      '+',
#'  "chr1", 150,    250,  3,      '-',
#'  "chr2", 1,      25,   4,      '+',
#'  "chr2", 200,    400,  5,      '-',
#'  "chr2", 400,    500,  6,      '+',
#'  "chr2", 450,    550,  7,      '+'
#' )
#' 
#' bed_merge(x)
#' 
#' bed_merge(x, max_dist = 100)
#' 
#' # merge intervals on same strand
#' bed_merge(dplyr::group_by(x, strand))
#' 
#' bed_merge(x, .value = sum(value))
#' 
#' @export
bed_merge <- function(x, max_dist = 0, ...) {
  
  x_groups <- groups(x)
  
  if (max_dist < 0)
    stop('max_dist must be positive', call. = FALSE)
  
  if ( ! is_sorted(x) ) {
    res <- bed_sort(x)
  } else {
    res <- x
  }
 
  res <- group_by(res, chrom, add = TRUE)
  
  res <- merge_impl(res, max_dist)
  
  dots <- list(.start = ~min(start), .end = ~max(end))
  dots <- c(dots, lazyeval::lazy_dots(...))
  
  res <- group_by_(res, .dots = c("chrom", ".id_merge", x_groups), add = TRUE)
    
  res <- summarize_(res, .dots = dots)
  res <- rename(res, start = .start, end = .end)
  res <- ungroup(res)
  # restore original grouping
  res <- group_by_(res, .dots = x_groups)
  res <- select(res, -.id_merge)
  res <- format_bed(res, x)
  
  attr(res, 'merged') <- TRUE

  res  
}

#' Ask whether a tbl is already merged.
#' 
#' @param x tbl of intervals
#' @noRd
is_merged <- function(x) {
  
  merged_attr <- attr(x, "merged")
  
  if ( is.null(merged_attr) || ! merged_attr ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
