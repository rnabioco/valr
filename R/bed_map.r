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
#' bed_map(bed_tbl, signal_tbl, sum = sum(value), max = max(value))
#' bed_map(bed_tbl, signal_tbl, concat = concat(value))
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
absmax <- function(.data) {
  abs(max(.data))
}

#' @export
#' @rdname bed_map
absmin <- function(.data) {
  abs(min(.data))
}

#' @export
#' @rdname bed_map
concat <- function(.data, sep = ',') {
  paste0(.data, collapse = sep)
}

#' @export
#' @rdname bed_map
distinct <- function(.data, sep = ',') {
  concat(unique(.data), sep = sep)
}

#' @export
#' @rdname bed_map
count <- function(.data) {
  length(.data)
}

#' @export
#' @rdname bed_map
count_distinct <- function(.data) {
  length(unique(.data))
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
