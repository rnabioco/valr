#' Subtract two sets of intervals.
#'
#' Subtract `y` intervals from `x` intervals.
#'
#' @param x [ivl_df]
#' @param y [ivl_df], or a path or URL to a bigWig (`.bw`) or bigBed (`.bb`)
#'   file. When a file is supplied, only the regions spanned by `x` are read
#'   from it (local files and `http(s)://` URLs are both supported), avoiding
#'   the cost of loading the entire file.
#' @param any remove any `x` intervals that overlap `y`
#'
#' @template min_overlap
#'
#' @template groups
#'
#' @family multiple set operations
#'
#' @seealso \url{https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 1,      100
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 50,     75
#' )
#'
#' bed_glyph(bed_subtract(x, y))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 100,    200,
#'   "chr1", 250,    400,
#'   "chr1", 500,    600,
#'   "chr1", 1000,   1200,
#'   "chr1", 1300,   1500
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 150,    175,
#'   "chr1", 510,    525,
#'   "chr1", 550,    575,
#'   "chr1", 900,    1050,
#'   "chr1", 1150,   1250,
#'   "chr1", 1299,   1501
#' )
#'
#' bed_subtract(x, y)
#'
#' bed_subtract(x, y, any = TRUE)
#'
#' @export
bed_subtract <- function(x, y, any = FALSE, min_overlap = 1L) {
  check_required(x)
  check_required(y)

  x <- check_interval(x)

  # `y` may be a path/URL to a bigWig/bigBed file; query only x's regions
  if (is_bigwig_path(y)) {
    y <- read_bigwig_regions(x, y)
  }

  y <- check_interval(y)

  # establish grouping with shared groups (and chrom)
  groups_xy <- shared_groups(x, y)
  groups_xy <- unique(c("chrom", groups_xy))

  # type convert grouping factors to characters if necessary and ungroup
  x <- convert_factors(x, groups_xy)
  y <- convert_factors(y, groups_xy)

  x <- group_by(x, across(all_of(groups_xy)))
  y <- group_by(y, across(all_of(groups_xy)))

  # find groups not in y
  not_y_grps <- setdiff(get_labels(x), get_labels(y))
  # keep x ivls from groups not found in y
  res_no_y <- semi_join(x, not_y_grps, by = colnames(not_y_grps))

  grp_indexes <- shared_group_indexes(x, y)

  if (any) {
    # collect and return x intervals without overlaps
    res <- intersect_impl(
      x,
      y,
      grp_indexes$x,
      grp_indexes$y,
      invert = TRUE,
      suffix_x = ".x",
      suffix_y = ".y",
      min_overlap = min_overlap
    )
    res <- tibble::as_tibble(res)
    anti <- filter(res, is.na(.data[[".overlap"]]))
    anti <- select(
      anti,
      all_of("chrom"),
      start = all_of("start.x"),
      end = all_of("end.x")
    )

    return(anti)
  }

  res <- subtract_impl(
    x,
    y,
    grp_indexes$x,
    grp_indexes$y,
    min_overlap
  )
  res <- tibble::as_tibble(res)
  res <- ungroup(res)
  res <- bind_rows(res, res_no_y)
  res <- bed_sort(res)

  res
}
