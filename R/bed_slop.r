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
bed_slop <- function(
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

  if (strand && !"strand" %in% colnames(x)) {
    cli::cli_abort("expected {.var strand} in {.var x}")
  }

  if (both != 0 && (left != 0 || right != 0)) {
    cli::cli_abort("ambiguous side spec for bed_slop")
  }

  if (fraction) {
    x <- mutate(x, .size = .data[["end"]] - .data[["start"]])
  }

  if (both != 0) {
    if (fraction) {
      res <- mutate(
        x,
        start = .data[["start"]] - round(both * .data[[".size"]]),
        end = .data[["end"]] + round(both * .data[[".size"]])
      )
    } else {
      res <- mutate(
        x,
        start = .data[["start"]] - both,
        end = .data[["end"]] + both
      )
    }
  } else {
    # calc left and right based on strand
    if (strand) {
      if (fraction) {
        res <- mutate(
          x,
          start = ifelse(
            .data[["strand"]] == "+",
            .data[["start"]] - round(left * .data[[".size"]]),
            .data[["start"]] - round(right * .data[[".size"]])
          ),
          end = ifelse(
            .data[["strand"]] == "+",
            .data[["end"]] + round(right * .data[[".size"]]),
            .data[["end"]] + round(left * .data[[".size"]])
          )
        )
      } else {
        res <- mutate(
          x,
          start = ifelse(
            .data[["strand"]] == "+",
            .data[["start"]] - left,
            .data[["start"]] - right
          ),
          end = ifelse(
            .data[["strand"]] == "+",
            .data[["end"]] + right,
            .data[["end"]] + left
          )
        )
      }
    } else {
      if (fraction) {
        res <- mutate(
          x,
          start = .data[["start"]] - round(left * .data[[".size"]]),
          end = .data[["end"]] + round(right * .data[[".size"]])
        )
      } else {
        res <- mutate(
          x,
          start = .data[["start"]] - left,
          end = .data[["end"]] + right
        )
      }
    }
  }

  if (fraction) {
    res <- select(res, -all_of(".size"))
  }

  res <- bound_intervals(res, genome, trim)

  res <- mutate(
    res,
    temp_start = .data[["start"]],
    temp_end = .data[["end"]]
  )

  res <- mutate(
    res,
    start = ifelse(
      .data[["temp_start"]] - .data[["temp_end"]] < 0,
      .data[["temp_start"]],
      .data[["temp_end"]]
    ),
    end = ifelse(
      .data[["temp_start"]] - .data[["temp_end"]] < 0,
      .data[["temp_end"]],
      .data[["temp_start"]]
    )
  )

  res <- select(res, -all_of(c("temp_start", "temp_end")))

  res
}
