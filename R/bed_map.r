#' Map signals over intervals.
#' 
#' Multiple mapping functions can be specified.
#' 
#' Variables from the `x` and `y` tables should be used with `.x` and `.y`
#' suffixes.
#' 
#' Existing grouping variables are maintained in the output, but have `.x` or 
#' `.y` suffixes`.
#' 
#' @param x tbl of intervals
#' @param y tbl of signals
#' @param ... name-value pairs of summary functions like \code{\link{min}()}, 
#'   \code{\link{count}()}, \code{\link{concat}()}.
#'   
#' @return \code{data_frame}
#'   
#' @seealso 
#' \url{http://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#' 
#' @note Column names have \code{.x} and \code{.y} suffixes.
#'   
#' @examples
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100, 250,
#'  "chr2", 250, 500)
#'  
#' y <- tibble::tribble(
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
#' bed_map(x, y, first(value.y))
#' bed_map(x, y, last(value.y))
#' 
#' bed_map(x, y, absmax = abs(max(value.y)))
#' bed_map(x, y, absmin = abs(min(value.y)))
#' bed_map(x, y, count = length(value.y))
#' bed_map(x, y, count_distinct = length(unique(value.y)))
#' 
#' bed_map(x, y, vals = values(value.y))
#' bed_map(x, y, vals.unique = values_unique(value.y))
#' 
#' @export
bed_map <- function(x, y, ...) {

  groups_x <- dplyr::groups(x)
  groups_y <- dplyr::groups(y)
  
  if('chrom' %in% c(groups_x, groups_y))
    stop('`chrom` cannot be used as grouping variable', call. = FALSE)
  
  if(!is.null(groups_x))
    groups_x <- stringr::str_c(groups_x, '.x')
  if(!is.null(groups_y))
    groups_y <- stringr::str_c(groups_y, '.y')
  
  res <- bed_intersect(x, y)
 
  groups_default <- c('chrom', 'start.x', 'end.x')
  res <- dplyr::group_by_(res, .dots = c(groups_default, groups_x, groups_y))
  
  res <- dplyr::summarize_(res, .dots = lazyeval::lazy_dots(...))
 
  # reassign original `x` groups 
  res <- dplyr::group_by_(res, .dots = c(groups_x))
 
  res 
}

#' @export
#' 
#' @param .data data
#' @param sep separator character
#' 
#' @rdname bed_map
concat <- function(.data, sep = ',') {
  paste0(.data, collapse = sep)
}

#' @export
#' @rdname bed_map
values_unique <- function(.data, sep = ',') {
  concat(unique(.data), sep = sep)
}

#' @export
#' @rdname bed_map
values <- function(.data, sep = ',') {
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
