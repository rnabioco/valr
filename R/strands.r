#' Flip strands in intervals.
#'
#' Flips positive (`+`) stranded intervals to negative (`-`) strands,
#' and vice-versa. Facilitates comparisons among intervals on opposing strands.
#'
#' @param x [ivl_df]
#'
#' @family utilities
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~strand,
#'   "chr1", 1,      100,  "+",
#'   "chr2", 1,      100,  "-"
#' )
#'
#' flip_strands(x)
#'
#' @export
flip_strands <- function(x) {
  if (!"strand" %in% colnames(x)) {
    cli::cli_abort("{.var strand} not found in {.var x}")
  }

  x <- check_interval(x)

  # remove existing groups
  groups_x <- groups(x)
  res <- ungroup(x)

  res <- mutate(res, .strand = ifelse(.data[["strand"]] == "+", "-", "+"))
  res <- select(res, -all_of("strand"))
  res <- select(res, everything(), strand = all_of(".strand"))

  res <- group_by(res, !!!groups_x)
  res
}
