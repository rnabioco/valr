#' Identify intersecting intervals.
#'
#' Report intersecting intervals from \code{x} and \code{y} tbls. Book-ended
#' intervals (or "touching" intervals) have \code{.overlap} values of \code{0})
#' in the output.
#'
#' @param x \code{\link{tbl_interval}}
#' @param ... single y  \code{\link{tbl_interval}}, multiple y \code{\link{tbl_interval}s}
#'                   or a list of y  \code{tbl_intervals}
#' @param invert report \code{x} intervals not in \code{y}
#' @param suffix colname suffixes in output
#'
#' @return \code{\link{tbl_interval}} with original columns from \code{x} and \code{y},
#'   suffixed with \code{.x} and \code{.y}, and a new \code{.overlap} column
#'   with the extent of overlap for the intersecting intervals.
#'
#'   If  multiple \code{y} tbls are supplied, then an additional column \code{source.y} will
#'   be reported that contains the variable names associated with each interval.
#'   If named tbls are supplied to \code{...} (i.e \code{a = y, b = z} or
#'    \code{list(a = y, b = z)}), then the supplied names will be reported instead
#'   of the variable names (see examples).
#'
#'
#' @template groups
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 25,      50,
#'   'chr1', 100,     125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75
#' )
#'
#' bed_glyph(bed_intersect(x, y))
#' bed_glyph(bed_intersect(x, y, invert = TRUE))
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    500,
#'   "chr2", 200,    400,
#'   "chr2", 300,    500,
#'   "chr2", 800,    900
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 150,    400,  100,
#'   "chr1", 500,    550,  100,
#'   "chr2", 230,    430,  200,
#'   "chr2", 350,    430,  300
#' )
#'
#' bed_intersect(x, y)
#'
#' bed_intersect(x, y, invert = TRUE)
#'
#' # start and end of each overlapping interval
#' res <- bed_intersect(x, y)
#' dplyr::mutate(res, start = pmax(start.x, start.y),
#'                    end = pmin(end.x, end.y))
#'
#' z <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 150,    400,  100,
#'   "chr1", 500,    550,  100,
#'   "chr2", 230,    430,  200,
#'   "chr2", 750,    900,  400
#' )
#'
#' bed_intersect(x, y, z)
#'
#' bed_intersect(x, exons = y, introns = z)
#'
#' # a list of tbl_intervals can also be passed
#' bed_intersect(x, list(exons = y, introns = z))
#'
#' @family multiple set operations
#' @seealso
#' \url{http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html}
#'
#'
#' @export
bed_intersect <- function(x, ..., invert = FALSE, suffix = c('.x', '.y')) {

  #determine if list supplied to ... or series of variables
  if (typeof(substitute(...)) == "symbol") {
    y_tbl <- list(...)
  } else {
    # extract out just a list not a list of lists
    y_tbl <- list(...)[[1]]
  }

  if (length(y_tbl) > 1){
    # multiple y tbls supplied
    if (is.null(names(y_tbl))){
      # name each tbl based on supplied variable
      dots <- eval(substitute(alist(...)))
      names(y_tbl) <- dots
    }
    #bind_rows preserves grouping
    y <- bind_rows(y_tbl, .id = "source")
    y <- select(y, -source, everything(), source)
  } else {
    # only one tbl supplied, so extract out single tbl from list
    y <- y_tbl[[1]]
  }

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_interval(y)) y <- tbl_interval(y)

  check_suffix(suffix)

  x <- group_by(x, chrom, add = TRUE)
  y <- group_by(y, chrom, add = TRUE)

  suffix <- list(x = suffix[1], y = suffix[2])

  res <- intersect_impl(x, y, suffix$x, suffix$y)

  if (invert) {
    colspec <- c('chrom' = 'chrom',
                 'start' = paste0('start', suffix$x),
                 'end' = paste0('end', suffix$x))
    res <- anti_join(x, res, by = colspec)
    res <- ungroup(res)
  }

  res
}
