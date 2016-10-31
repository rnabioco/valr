#' Map signals over intervals.
#' 
#' Multiple mapping functions can be specified.
#' 
#' @inheritParams bed_intersect
#' @param min_overlap minimum overlap for intervals.
#' @param ... name-value pairs of summary functions like \code{\link{min}()}, 
#'   \code{\link{count}()}, \code{\link{concat}()}.
#'   
#' @return \code{data_frame}
#'   
#' @family multi-set-ops
#' @seealso 
#' \url{http://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#' 
#' @note Book-ended intervals are not reported by default, but can be included
#'   by setting \code{min_overlap} to \code{0}.
#'   
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1',      1,      100
#' )
#' 
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1',      1,     20,    10,
#'   'chr1',      30,    50,    20,
#'   'chr1',      90,    120,   30
#' )
#' 
#' bed_glyph(bed_map(x, y, value = sum(value)), label = 'value')
#' 
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
#' # also mean, median, sd etc
#' bed_map(x, y, .sum = sum(value))
#' bed_map(x, y, .min = min(value), .max = max(value))
#' 
#' bed_map(x, y, .concat = concat(value))
#' 
#' # can also use `nth` family from dplyr
#' bed_map(x, y, .first = dplyr::first(value))
#' bed_map(x, y, .last = dplyr::last(value))
#' 
#' bed_map(x, y, .absmax = abs(max(value)))
#' bed_map(x, y, .absmin = abs(min(value)))
#' bed_map(x, y, .count = length(value))
#' bed_map(x, y, .count_distinct = length(unique(value)))
#' 
#' bed_map(x, y, .vals = values(value))
#' bed_map(x, y, .vals.unique = values_unique(value))
#' 
#' @export
bed_map <- function(x, y, ..., invert = FALSE,
                    suffix = c('.x', '.y'),
                    min_overlap = 1) {
  
  groups_x <- groups(x)
  groups_y <- groups(y)
  
  # used only to get the `x` suffix; `y` suffix is ignored` 
  suffix <- list(x = suffix[1], y = suffix[2])
 
  # `x` names are suffixed to use for grouping later 
  x_names <- colnames(x)[!colnames(x) %in% "chrom"]
  x_names_suffix <- stringr::str_c(x_names, suffix$x)
 
  # note that `y` columns have no suffix so can be referred to by the original names
  res <- bed_intersect(x, y, invert = invert, suffix = c(suffix$x, ''))
  
  res <- filter(res, .overlap >= min_overlap)
  
  ##  map supplied functions to each set of intervals
  res <- group_by_(res, .dots = c("chrom", x_names_suffix))
  res <- summarize_(res, .dots = lazyeval::lazy_dots(...))
  res <- ungroup(res)
  ## remove x suffix, but don't pattern match with '.' regex
  names_no_x <- stringr::str_replace(names(res), stringr::fixed(suffix$x), '')
  names(res) <- names_no_x
  
  # find rows of `x` that did not intersect
  x_not <- anti_join(x, res, by = c("chrom", x_names))
  
  res <- bind_rows(res, x_not)
  res <- bed_sort(res)
 
  # reassign original `x` groups. `y` groups are gone at this point
  res <- group_by_(res, .dots = c(groups_x))
  
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
