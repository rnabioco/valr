#' Create flanking intervals from input intervals.
#'
#' @param x [tbl_interval()]
#' @param genome [tbl_genome()]
#' @param both number of bases on both sizes
#' @param left number of bases on left side
#' @param right number of bases on right side
#' @param fraction define flanks based on fraction of interval length
#' @param strand define `left` and `right` based on strand
#' @param trim adjust coordinates for out-of-bounds intervals
#' @param ... extra arguments (not used)
#'
#' @return [tbl_interval()]
#'
#' @family single set operations
#' @seealso
#'   \url{http://bedtools.readthedocs.org/en/latest/content/tools/flank.html}
#'
#' @examples
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1',      25,      50,
#'   'chr1',      100,     125
#' )
#'
#' genome <- trbl_genome(
#'   ~chrom, ~size,
#'   'chr1', 130
#' )
#'
#' bed_glyph(bed_flank(x, genome, both = 20))
#'
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'  "chr1", 500,    1000, '.',   '.',    '+',
#'  "chr1", 1000,   1500, '.',   '.',    '-'
#' )
#'
#' genome <- trbl_genome(
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

bed_flank <- function(x, genome, both = 0, left = 0,
                      right = 0, fraction = FALSE,
                      strand = FALSE, trim = FALSE, ...) {

  if (!is.tbl_interval(x)) x <- tbl_interval(x)

  res <- flank_impl(x, genome, both, left,
                    right, fraction, strand, trim)

  arrange(res, chrom, start)

}

