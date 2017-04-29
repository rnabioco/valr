#' Provide working directory for valr example files.
#'
#' @param path path to file
#'
#' @examples
#' valr_example('hg19.chrom.sizes.gz')
#'
#' @export
valr_example <- function(path) {
  # https://twitter.com/JennyBryan/status/780150538654527488
  system.file("extdata", path, package = 'valr', mustWork = TRUE)
}

#' reformat tbl column ordering based upon another tbl
#'
#' `reorder_names` returns a tbl whose columns are ordered by another tbl.
#' The `x` tbl columns are reordered based on the `y` columns
#' ordering. `x` columns that do not exist in `y` are moved to the
#' last column.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#'
#' @examples
#' # names out of order
#' x <- trbl_interval(
#'   ~end, ~chrom, ~start, ~value,
#'   75,   "chr1", 125,    10
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end,  ~scores,
#'   "chr1", 50,     100,   1.2,
#'   "chr1", 100,    150,   2.4
#' )
#'
#' reorder_names(x, y)
#' @noRd
reorder_names <- function(x, y) {
  names_x <- names(x)
  names_y <- names(y)

  names_x <- names_x[order(match(names_x,names_y))]

  x <- select(x, one_of(names_x))
  x
}


#' Identify groups shared between to tbls
#'
#' Identify minimum shared groups between `x` and `y` tbls. Returns
#' `NULL` if there are no shared groups.
#'
#' @param x [tbl_interval()]
#' @param y [tbl_interval()]
#'
#' @return `list` of groups or `NULL`
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 150,    400,  100,
#'   "chr2", 230,    430,  200
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 50,     100,  1,
#'   "chr1", 100,    150,  2
#' )
#'
#' x <- dplyr::group_by(x, chrom, value)
#' y <- dplyr::group_by(y, chrom, value)
#' shared_groups(x, y)
#'
#' y <- dplyr::group_by(y, chrom)
#' shared_groups(x, y)
#'
#' y <- dplyr::ungroup(y)
#' shared_groups(x, y)
#' @noRd
shared_groups <- function(x, y) {

  groups_x <- groups(x)
  groups_y <- groups(y)

  groups_xy <- intersect(groups_x, groups_y)
  if (length(groups_xy) == 0){
    groups_xy <- NULL
  }
  groups_xy
}

# dplyr::check_suffix
check_suffix <- function(suffix) {
  if (!is.character(suffix) || length(suffix) != 2)
    stop("`suffix` must be a character vector of length 2.", call. = FALSE)
}

#' Return group labels from tbl_df
#' @param grp_df grouped tbl_df
#' @return `tibble` of grouping labels or `NULL` if no groups present
#' @noRd
get_labels <- function(grp_tbl) attr(grp_tbl, "labels")
