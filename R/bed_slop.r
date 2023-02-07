#' Increase the size of input intervals.
#'
#' @inheritParams bed_flank
#'
#' @return [ivl_df]
#'
#' @family single set operations
#'
#' @seealso
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/slop.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 110,    120,
#'   "chr1", 225,    235
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 400
#' )
#'
#' bed_glyph(bed_slop(x, genome, both = 20, trim = TRUE))
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 5000
#' )
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~name, ~score, ~strand,
#'   "chr1", 500, 1000, ".", ".", "+",
#'   "chr1", 1000, 1500, ".", ".", "-"
#' )
#'
#' bed_slop(x, genome, left = 100)
#'
#' bed_slop(x, genome, right = 100)
#'
#' bed_slop(x, genome, both = 100)
#'
#' bed_slop(x, genome, both = 0.5, fraction = TRUE)
#'
#' @export
bed_slop <- function(x, genome, both = 0, left = 0,
                     right = 0, fraction = FALSE,
                     strand = FALSE, trim = FALSE, ...) {
  check_required(x)
  check_required(genome)

  x <- check_interval(x)
  genome <- check_genome(genome)

  if (strand && !"strand" %in% colnames(x)) {
    cli::cli_abort("expected {.var strand} in {.var x}")
  }

  if (both != 0 && (left != 0 || right != 0)) {
    cli::cli_abort("ambiguous side spec for bed_slop")
  }

  if (fraction) x <- mutate(x, .size = end - start)

  if (both != 0) {
    if (fraction) {
      res <- mutate(
        x,
        start = start - round(both * .size),
        end = end + round(both * .size)
      )
    } else {
      res <- mutate(
        x,
        start = start - both,
        end = end + both
      )
    }
  } else {
    # calc left and right based on strand
    if (strand) {
      if (fraction) {
        res <- mutate(
          x,
          start = ifelse(strand == "+",
            start - round(left * .size),
            start - round(right * .size)
          ),
          end = ifelse(strand == "+",
            end + round(right * .size),
            end + round(left * .size)
          )
        )
      } else {
        res <- mutate(
          x,
          start = ifelse(strand == "+",
            start - left,
            start - right
          ),
          end = ifelse(strand == "+",
            end + right,
            end + left
          )
        )
      }
    } else {
      if (fraction) {
        res <- mutate(
          x,
          start = start - round(left * .size),
          end = end + round(right * .size)
        )
      } else {
        res <- mutate(
          x,
          start = start - left,
          end = end + right
        )
      }
    }
  }

  if (fraction) res <- select(res, -.size)

  res <- bound_intervals(res, genome, trim)
  res <- bed_sort(res)

  res <- mutate(res,
    temp_start = start,
    temp_end = end
  )

  res <- mutate(res,
    start = ifelse(temp_start - temp_end < 0,
      temp_start, temp_end
    ),
    end = ifelse(temp_start - temp_end < 0,
      temp_end, temp_start
    )
  )

  res <- select(res, -temp_start, -temp_end)

  res
}
