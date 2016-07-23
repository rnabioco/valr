#' Map signals over intervals.
#' 
#' @param x tbl of intervals
#' @param y tbl of signals
#' @param add_group additional grouping vars
#' @param ... name-value pairs of summary functions like \code{\link{min}()}, 
#'   \code{\link{count}()}, \code{\link{concat}()}. 
#'   
#' @return \code{data_frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#' 
#' @note Column names have \code{.x} and \code{.y} suffixes.
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
#' bed_map(x, y, sum = sum(value.y))
#' bed_map(x, y, min = min(value.y), max = max(value.y))
#' 
#' bed_map(x, y, concat(value.y))
#' bed_map(x, y, distinct(value.y))
#' bed_map(x, y, first(value.y))
#' bed_map(x, y, last(value.y))
#' 
#' bed_map(x, y, absmax = abs(max(value.y)))
#' bed_map(x, y, absmin = abs(min(value.y)))
#' bed_map(x, y, count = length(value.y))
#' bed_map(x, y, count_distinct = length(unique(value.y)))
#' 
#' # use decreasing = TRUE to reverse numbers
#' bed_map(x, y, distinct_num = distinct(sort(value.y)))
#' 
#' @export
bed_map <- function(x, y, ..., add_group = NULL) {

  groups_default <- c('chrom', 'start.x', 'end.x')
  res <- bed_intersect(x, y)
  res <- group_by_(res, .dots = c(groups_default, add_group))
  res <- summarize_(res, .dots = lazyeval::lazy_dots(...))

  res 
}

#' @export
#' @rdname bed_map
#' 
#' @param .data vector of things
#' @param sep separator character
#' 
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
