#' Sort a set of intervals.
#'
#' @param x [ivl_df]
#' @param by_size sort by interval size
#' @param by_chrom sort within chromosome
#' @param reverse reverse sort order
#'
#' @seealso
#' \url{https://bedtools.readthedocs.io/en/latest/content/tools/sort.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr8", 500,    1000,
#'   "chr8", 1000,   5000,
#'   "chr8", 100,    200,
#'   "chr1", 100,    300,
#'   "chr1", 100,    200
#' )
#'
#' # sort by chrom and start
#' bed_sort(x)
#'
#' # reverse sort order
#' bed_sort(x, reverse = TRUE)
#'
#' # sort by interval size
#' bed_sort(x, by_size = TRUE)
#'
#' # sort by decreasing interval size
#' bed_sort(x, by_size = TRUE, reverse = TRUE)
#'
#' # sort by interval size within chrom
#' bed_sort(x, by_size = TRUE, by_chrom = TRUE)
#'
#' @export
bed_sort <- function(x, by_size = FALSE, by_chrom = FALSE, reverse = FALSE) {
  check_required(x)
  x <- check_interval(x)

  if (by_size) {
    res <- mutate(x, .size = .data[["end"]] - .data[["start"]])

    if (by_chrom) {
      if (reverse) {
        res <- res[order(res$chrom, -res$.size, method = "radix"), ]
      } else {
        res <- res[order(res$chrom, res$.size, method = "radix"), ]
      }
    } else {
      if (reverse) {
        res <- res[order(-res$.size, method = "radix"), ]
      } else {
        res <- res[order(res$.size, method = "radix"), ]
      }
    }

    # remove .size column and groups in result
    res <- select(res, -all_of(".size"))
  } else {
    # sort by coordinate
    if (reverse) {
      res <- x[order(x$chrom, -x[["start"]], method = "radix"), ]
    } else {
      res <- x[order(x$chrom, x[["start"]], x[["end"]], method = "radix"), ]
    }
  }

  res
}
