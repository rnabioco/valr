#' Adjust intervals by a fixed size.
#'
#' Out-of-bounds intervals are removed by default.
#'
#' @param x [ivl_df]
#' @param genome [ivl_df]
#' @param size number of bases to shift. positive numbers shift right, negative shift left.
#' @param fraction define `size` as a fraction of interval
#' @param trim adjust coordinates for out-of-bounds intervals
#'
#' @return [ivl_df]
#'
#' @family single set operations
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/shift.html}
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
#'   "chr1", 125
#' )
#'
#' bed_glyph(bed_shift(x, genome, size = -20))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~strand,
#'   "chr1", 100,    150,  "+",
#'   "chr1", 200,    250,  "+",
#'   "chr2", 300,    350,  "+",
#'   "chr2", 400,    450,  "-",
#'   "chr3", 500,    550,  "-",
#'   "chr3", 600,    650,  "-"
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 1000,
#'   "chr2", 2000,
#'   "chr3", 3000
#' )
#'
#' bed_shift(x, genome, 100)
#'
#' bed_shift(x, genome, fraction = 0.5)
#'
#' # shift with respect to strand
#' stranded <- dplyr::group_by(x, strand)
#' bed_shift(stranded, genome, 100)
#'
#' @export
bed_shift <- function(x, genome, size = 0, fraction = 0, trim = FALSE) {
  check_required(x)
  check_required(genome)

  x <- check_interval(x)
  genome <- check_genome(genome)

  stranded <- "strand" %in% groups(x)

  # shift invervals
  if (!stranded && !fraction) {
    res <- mutate(
      x,
      start = .data[["start"]] + size,
      end = .data[["end"]] + size
    )
  }

  # shift by percent of interval size
  if (!stranded && fraction) {
    res <- mutate(x, .size = .data[["end"]] - .data[["start"]])
    res <- mutate(
      res,
      start = .data[["start"]] + round(.data[[".size"]] * fraction),
      end = .data[["end"]] + round(.data[[".size"]] * fraction)
    )
    res <- select(res, -all_of(".size"))
  }

  # shift by strand
  if (stranded && !fraction) {
    res <- mutate(
      x,
      start = ifelse(
        .data[["strand"]] == "+",
        .data[["start"]] + size,
        .data[["start"]] - size
      ),
      end = ifelse(
        .data[["strand"]] == "+",
        .data[["end"]] + size,
        .data[["end"]] - size
      )
    )
  }

  # shift by strand and percent
  if (stranded && fraction) {
    res <- mutate(x, .size = .data[["end"]] - .data[["start"]])
    res <- mutate(
      res,
      start = ifelse(
        .data[["strand"]] == "+",
        .data[["start"]] + round(.data[[".size"]] * fraction),
        .data[["start"]] - round(.data[[".size"]] * fraction)
      ),
      end = ifelse(
        .data[["strand"]] == "+",
        .data[["end"]] + round(.data[[".size"]] * fraction),
        .data[["end"]] - round(.data[[".size"]] * fraction)
      )
    )
    res <- select(res, -all_of(".size"))
  }

  res <- bound_intervals(res, genome, trim)

  res
}
