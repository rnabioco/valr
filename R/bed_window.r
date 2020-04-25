#' Identify intervals within a specified distance.
#'
#' @param x [ivl_df]
#' @param y  [ivl_df]
#' @param ... params for bed_slop and bed_intersect
#' @inheritParams bed_slop
#' @inheritParams bed_intersect
#'
#' @template groups
#'
#' @family multiple set operations
#' @examples
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  'chr1', 25,     50,
#'  'chr1', 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 60,     75
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   'chr1', 125
#' )
#'
#' bed_glyph(bed_window(x, y, genome, both = 15))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 10,    100,
#'   "chr2", 200,    400,
#'   "chr2", 300,    500,
#'   "chr2", 800,    900
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 150,    400,
#'   "chr2", 230,    430,
#'   "chr2", 350,    430
#' )
#'
#' genome <- tibble::tribble(
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
  x <- check_interval(x)
  y <- check_interval(y)
  genome <- check_genome(genome)

  x <- mutate(x, .start = start, .end = end)

  # capture command line args
  cmd_args <- list(...)

  # get arguments for bed_slop and bed_intersect
  slop_arg_names <- names(formals(bed_slop))
  intersect_arg_names <- names(formals(bed_intersect))

  # parse supplied args into those for bed_slop or bed_intersect
  slop_args <- cmd_args[names(cmd_args) %in% slop_arg_names]
  intersect_args <- cmd_args[names(cmd_args) %in% intersect_arg_names]

  # pass new list of args to bed_slop
  slop_x <- do.call(
    bed_slop,
    c(list("x" = x, "genome" = genome), slop_args)
  )

  # pass new list of args to bed_intersect
  res <- do.call(
    bed_intersect,
    c(list("x" = slop_x, "y" = y), intersect_args)
  )

  res <- mutate(res, start.x = .start.x, end.x = .end.x)

  res <- ungroup(res)

  res <- select(res, -.start.x, -.end.x)

  res
}
