#' Create flanking intervals from input intervals.
#'
#' @param x [ivl_df]
#' @param genome [genome_df]
#' @param both number of bases on both sizes
#' @param left number of bases on left side
#' @param right number of bases on right side
#' @param fraction define flanks based on fraction of interval length
#' @param strand define `left` and `right` based on strand
#' @param trim adjust coordinates for out-of-bounds intervals
#' @param ... extra arguments (not used)
#'
#' @return [ivl_df]
#'
#' @family single set operations
#' @seealso
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/flank.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 25, 50,
#'   "chr1", 100, 125
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 130
#' )
#'
#' bed_glyph(bed_flank(x, genome, both = 20))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 500,    1000, ".",   ".",    "+",
#'   "chr1", 1000,   1500, ".",   ".",    "-"
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 5000
#' )
#'
#' bed_flank(x, genome, left = 100)
#'
#' bed_flank(x, genome, right = 100)
#'
#' bed_flank(x, genome, both = 100)
#'
#' bed_flank(x, genome, both = 0.5, fraction = TRUE)
#'
#' @export

bed_flank <- function(
  x,
  genome,
  both = 0,
  left = 0,
  right = 0,
  fraction = FALSE,
  strand = FALSE,
  trim = FALSE,
  ...
) {
  check_required(x)
  check_required(genome)

  x <- check_interval(x)
  genome <- check_genome(genome)

  if (!any(c(both, left, right) > 0)) {
    cli::cli_abort(
      "one of {.var both}, {.var left}, or {.var right} must be a positive value"
    )
  }

  if (strand && !"strand" %in% colnames(x)) {
    cli::cli_abort("expected {.var strand} column in {.var x}")
  }

  if (both != 0 && (left != 0 || right != 0)) {
    cli::cli_abort("ambiguous side spec for bed_flank")
  }

  if (both) {
    left <- right <- both
  }

  res <- flank_impl(
    x,
    genome,
    both,
    left,
    right,
    fraction,
    strand,
    trim
  )

  res <- bed_sort(res)

  res
}
