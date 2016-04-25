#' Merge overlapping intervals.
#'
#' @param bed_tbl tbl of intervals 
#' @param max_dist maximum distance between intervals to merge
#' @param .keep keep dot columns (\code{.merge_id}, \code{.overlap})
#' @param ... name-value pairs that specify merging operations
#'
#' @return \code{data_frame}
#'  
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/merge.html}
#' 
#' @examples 
#' bed_tbl <- tibble::frame_data(
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
#' bed_merge(bed_tbl)
#' bed_merge(bed_tbl, max_dist = 100)
#' bed_merge(bed_tbl, strand = TRUE)
#' bed_merge(bed_tbl, .value = sum(value))
#' 
#' @export
bed_merge <- function(bed_tbl, max_dist = 0, strand = FALSE, ...) {
 
  assert_that(max_dist >= 0)
  
  if ( ! is_sorted(bed_tbl) ) {
    res <- bed_sort(bed_tbl)
  }
 
  res <- group_by(res, chrom)
  
  if (strand)
    res <- group_by(res, strand, add = TRUE)

  res <- merge_impl(res, max_dist)
  
  dots <- list(.start = ~min(start), .end = ~max(end))
  dots <- c(dots, lazyeval::lazy_dots(...))
  
  res <- res %>% 
    group_by(chrom, .merge_id) %>%
    summarize_(.dots = dots) %>%
    rename(start = .start, end = .end) %>%
    select(-.merge_id)
  
  attr(res, 'merged') <- TRUE

  res  
}

#' Ask whether a tbl is already merged.
#' 
#' @export
is_merged <- function(bed_tbl) {
  
  merged_attr <- attr(bed_tbl, "merged")
  
  if ( is.null(merged_attr) || ! merged_attr ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

