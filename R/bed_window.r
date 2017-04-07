#' Identify intervals within a specified distance.
#'
#' @param ... params for bed_slop and bed_intersect
#' @inheritParams bed_slop
#' @inheritParams bed_intersect
#'
#' @template groups
#'
#' @family multiple set operations
#' @examples
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  'chr1', 25,     50,
#'  'chr1', 100,    125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 60,     75
#' )
#'
#' genome <- trbl_genome(
#'   ~chrom, ~size,
#'   'chr1', 125
#' )
#'
#' bed_glyph(bed_window(x, y, genome, both = 15))
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10,    100,
#'   "chr2", 200,    400,
#'   "chr2", 300,    500,
#'   "chr2", 800,    900
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   "chr1", 150,    400,
#'   "chr2", 230,    430,
#'   "chr2", 350,    430
#' )
#'
#' genome <- trbl_interval(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#'
#' bed_window(x, y, genome, both = 100)
#'
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/window.html}
#'
#' @export
bed_window <- function(x, y, genome, ...) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)
  if (!is.tbl_interval(y)) y <- tbl_interval(y)
  if (!is.tbl_genome(genome)) genome <- tbl_genome(genome)

  x <- mutate(x, .start = start, .end = end)

  slop_x <- bed_slop(x, genome, ...)

  res <- bed_intersect(slop_x, y, ...)

  res <- mutate(res, start.x = .start.x, end.x = .end.x)

  res <- ungroup(res)

  res <- select(res, -.start.x, -.end.x)

  res
}
