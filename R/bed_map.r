#' Map signals over intervals
#' 
#' @param bed_tbl tbl of intervals 
#' @param signal_tbl tbl of signals 
#' @param ... name-value pairs of summary functions like \code{\link{min}()},
#'   \code{\link{count}()}, \code{\link{concat}()}
#'        
#' @return \code{data_frame}
#' 
#' @examples
#' bed_tbl <- tibble::frame_data(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100, 250,
#'  "chr2", 250, 500)
#'  
#' signal_tbl <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~value,
#'  "chr1", 100, 250, 10,
#'  "chr1", 150, 250, 20,
#'  "chr2", 250, 500, 500)
#' 
#' # mean, median, sd etc
#' bed_map(bed_tbl, signal_tbl, sum = sum(value))
#' bed_map(bed_tbl, signal_tbl, min = min(value), max = max(value))
#' 
#' bed_map(bed_tbl, signal_tbl, concat(value))
#' bed_map(bed_tbl, signal_tbl, distinct(value))
#' bed_map(bed_tbl, signal_tbl, first(value))
#' bed_map(bed_tbl, signal_tbl, last(value))
#' 
#' bed_map(bed_tbl, signal_tbl, absmax = abs(max(value)))
#' bed_map(bed_tbl, signal_tbl, absmin = abs(min(value)))
#' bed_map(bed_tbl, signal_tbl, count = length(value))
#' bed_map(bed_tbl, signal_tbl, count_distinct = length(unique(value)))
#' 
#' # use decreasing = TRUE to reverse numbers
#' bed_map(bed_tbl, signal_tbl, distinct_num = distinct(sort(value)))
#' 
#' @export
bed_map <- function(bed_tbl, signal_tbl, ...) {
 
  res <- bed_intersect(bed_tbl, signal_tbl) %>%
    group_by(chrom, start.x, end.x) %>%
    summarize_(.dots = lazyeval::lazy_dots(...)) %>%
    rename(start = start.x, end = end.x) %>%
    ungroup()

  res 
}

#' @export
#' @rdname bed_map
concat <- function(.data, sep = ',') {
  paste0(.data, collapse = sep)
}

#' @export
#' @rdname bed_map
distinct_only <- function(.data, sep = ',') {
  concat(unique(.data), sep = sep)
}

#' @export
#' @rdname bed_map
distinct <- function(.data, sep = ',') {
  concat(rle(.data)$values, sep = sep)
}

#' @export
#' @rdname bed_map
first <- function(.data) {
  head(.data, n = 1)
}

#' @export
#' @rdname bed_map
last <- function(.data) {
  tail(.data, n = 1)
}
