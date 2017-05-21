#' Flip strands in intervals.
#'
#' Flips positive (`+`) stranded intervals to negative (`-`) strands,
#' and vice-versa. Facilitates comparisons among intervals on opposing strands.
#'
#' @param x [tbl_interval()]
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~strand,
#'   'chr1', 1,      100,  '+',
#'   'chr2', 1,      100,  '-'
#' )
#'
#' flip_strands(x)
#'
#' @export
flip_strands <- function(x) {

  if (! 'strand' %in% colnames(x))
    stop('`strand` column not found in `x`', call. = FALSE)

  if (!is.tbl_interval(x)) x <- tbl_interval(x)

  # remove existing groups
  groups_x <- groups(x)
  res <- ungroup(x)

  res <- mutate(res, .strand = ifelse(strand == '+', '-', '+'))
  res <- select(res, -strand)
  res <- rename(res, strand = .strand)

  res <- group_by(res, !!! groups_x)
  res
}
