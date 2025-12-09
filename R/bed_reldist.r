#' Compute relative distances between intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param detail report relative distances for each `x` interval.
#'
#' @family interval statistics
#'
#' @return
#' If `detail = FALSE`, a [ivl_df] that summarizes
#' calculated `.reldist` values with the following columns:
#'
#'   - `.reldist` relative distance metric
#'   - `.counts` number of metric observations
#'   - `.total` total observations
#'   - `.freq` frequency of observation
#'
#' If `detail = TRUE`, the `.reldist` column reports the relative
#' distance for each input `x` interval.
#'
#' @template stats
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/reldist.html}
#'
#' @examples
#' genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
#'
#' x <- bed_random(genome, seed = 1010486)
#' y <- bed_random(genome, seed = 9203911)
#'
#' bed_reldist(x, y)
#'
#' bed_reldist(x, y, detail = TRUE)
#'
#' @export
bed_reldist <- function(x, y, detail = FALSE) {
  check_required(x)
  check_required(y)

  x <- check_interval(x)
  y <- check_interval(y)

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(as.character(c("chrom", groups_xy)))
  groups_vars <- rlang::syms(groups_xy)

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, !!!groups_vars)
  y <- group_by(y, !!!groups_vars)

  grp_indexes <- shared_group_indexes(x, y)

  res <- dist_impl(x, y, grp_indexes$x, grp_indexes$y, distcalc = "reldist")
  res <- tibble::as_tibble(res)

  if (detail) {
    return(res)
  }

  res[[".reldist"]] <- floor(res[[".reldist"]] * 100) / 100
  nr <- nrow(res)
  res <- group_by(res, .data[[".reldist"]])
  res <- summarize(
    res,
    .counts = n(),
    .total = nr,
    .freq = .data[[".counts"]] / .data[[".total"]]
  )
  res
}
