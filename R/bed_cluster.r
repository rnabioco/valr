#' Cluster neighboring intervals.
#'
#' The ouput `.id` column can be used in downstream grouping operations. Default
#' `max_dist = 0` means that both overlapping and book-ended intervals will be
#' clustered.
#'
#' @param x [ivl_df]
#' @param max_dist maximum distance between clustered intervals.
#'
#' @template groups
#'
#' @return [ivl_df] with `.id` column specifying sets of clustered intervals.
#'
#' @family single set operations
#'
#' @seealso
#' \url{http://bedtools.readthedocs.org/en/latest/content/tools/cluster.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 100,    200,
#'   'chr1', 180,    250,
#'   'chr1', 250,    500,
#'   'chr1', 501,    1000,
#'   'chr2', 1,      100,
#'   'chr2', 150,    200
#' )
#'
#' bed_cluster(x)
#'
#' # glyph illustrating clustering of overlapping and book-ended intervals
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 1,      10,
#'   'chr1', 5,      20,
#'   'chr1', 30,     40,
#'   'chr1', 40,     50,
#'   'chr1', 80,     90
#' )
#'
#' bed_glyph(bed_cluster(x), label = '.id')
#'
#' @export
bed_cluster <- function(x, max_dist = 0) {
  x <- check_interval(x)

  groups <- rlang::syms(unique(c("chrom", group_vars(x))))
  res <- group_by(x, !!! groups)

  res <- bed_sort(res)

  res <- merge_impl(res, max_dist, collapse = FALSE)

  res <- mutate(res, .id = .id_merge)
  res <- select(res, !! quo(-one_of(".id_merge", ".overlap_merge")))
  res <- ungroup(res)
  res <- mutate(
    res,
    .id = match(.id, unique(.id))
  )
  res
}
