#' Map signals over intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of signals
#' @param ... name-value pairs of summary functions like \code{\link{min}()}, 
#'   \code{\link{count}()}, \code{\link{concat}()}. \code{start} and \code{end}
#'   colnames in returned data from .x suffixes.
#'   
#' @return \code{data_frame}
#'   
#' @examples
#' x <- tibble::frame_data(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100, 250,
#'  "chr2", 250, 500)
#'  
#' y <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~value,
#'  "chr1", 100, 250, 10,
#'  "chr1", 150, 250, 20,
#'  "chr2", 250, 500, 500)
#' 
#' # mean, median, sd etc
#' bed_map(x, y, sum = sum(value))
#' bed_map(x, y, min = min(value), max = max(value))
#' 
#' bed_map(x, y, concat(value))
#' bed_map(x, y, distinct(value))
#' bed_map(x, y, first(value))
#' bed_map(x, y, last(value))
#' 
#' bed_map(x, y, absmax = abs(max(value)))
#' bed_map(x, y, absmin = abs(min(value)))
#' bed_map(x, y, count = length(value))
#' bed_map(x, y, count_distinct = length(unique(value)))
#' 
#' # use decreasing = TRUE to reverse numbers
#' bed_map(x, y, distinct_num = distinct(sort(value)))
#' 
#' @export
bed_map <- function(x, y, ...) {

  res <- bed_intersect(x, y) %>%
    group_by(chrom, start.x, end.x) %>%
    summarize_(.dots = lazyeval::lazy_dots(...))

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
