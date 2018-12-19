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
  system.file("extdata", path, package = "valr", mustWork = TRUE)
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

  names_x <- names_x[order(match(names_x, names_y))]

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
  if (length(groups_xy) == 0) {
    groups_xy <- NULL
  }
  groups_xy
}

# dplyr::check_suffix
check_suffix <- function(suffix) {
  if (!is.character(suffix) || length(suffix) != 2) {
    stop("`suffix` must be a character vector of length 2.", call. = FALSE)
  }
}

#' Return group labels from tbl_df
#' @param grp_df grouped tbl_df
#' @return `tibble` of grouping labels or `NULL` if no groups present
#' @noRd
get_labels <- function(grp_tbl) {

  dplyr_version <- utils::packageVersion("dplyr")

  if (dplyr_version < "0.7.9.9000") {
    grp_df <- attr(grp_tbl, "labels")
  } else if (dplyr_version >= "0.7.9.9000") {
    grp_df <- attr(grp_tbl, "groups")
    grp_df <- grp_df[, !colnames(grp_df) %in% ".rows"]
  } else {
    stop("unable to extract labels from input dataframe")
  }
  grp_df
}

#' Update grouped dataframe attributes for dplyr 0.8.0 version
#' @param df grouped data frame
#' @return grouped data frame with attributes matching dplyr 0.8.0 version
#' @noRd
update_groups <- function(df){
  r_indexes <- lapply(attr(df, "indices"),
                      function(x) x + 1L)

  attr(df, "groups") <- as_tibble(attr(df, "labels"))
  attr(df, "groups")$.rows <- r_indexes
  attr(df, "biggest_group_size") <- NULL
  attr(df, "group_sizes") <- NULL
  attr(df, "indices") <- NULL
  attr(df, "labels") <- NULL
  attr(df, "vars") <- NULL
  attr(df, "drop") <- NULL
  df
}

#' Type convert factors if they are grouping columns
#' @param x data frame
#' @param group_cols group columns to type convert if factors
#' @return ungrouped dataframe
#' @noRd
convert_factors <- function(x, group_cols){

  contains_factor <- sapply(x[, group_cols], is.factor)
  if (any(contains_factor)){

    warning("Factors are not allowed for grouping: ",
            paste(group_cols[contains_factor],
                  collapse = ","),
            " will be treated as characters",
            call. = FALSE)

    x <- ungroup(x)
    convert_cols <- group_cols[contains_factor]
    x <- mutate_at(x, .vars = convert_cols, as.character)
  } else {
    ungroup(x)
  }
  x
}


#' Get indexes of groups in each grouped data.frame that are found in the other data.frame
#'
#' @param x grouped data.frame
#' @param y grouped data.frame
#' @return named list with integer vector of indexes of groups shared between data.frames
#' @noRd
shared_group_indexes <- function(x, y){
  x <- get_group_data(x)
  y <- get_group_data(y)
  shared_rows(x, y)
}

#' Get indexes of rows in each data.frame that are found in the other data.frame
#'
#' By default only columns with the same names in the two data.frames will
#' be compared
#'
#' @param x data.frame
#' @param y data.frame
#' @return named list with integer vector of indexes shared between data.frames
#' @noRd
shared_rows <- function(x, y){
  # based on plyr::match_df
  shared_cols <- intersect(colnames(x), colnames(y))
  combined_df <- bind_rows(x[shared_cols], y[shared_cols])

  keys <- id(combined_df, drop = TRUE)
  n_x <- nrow(x)
  n_y <- nrow(y)
  keys <- list(x = keys[seq_len(n_x)],
               y = keys[n_x + seq_len(n_y)],
               n = attr(keys, "n"))
  x_indexes <- which(keys$x %in% keys$y)
  y_indexes <- which(keys$y %in% keys$x)
  list(x = x_indexes,
       y = y_indexes)
}

#' Get data.frame from groups attribute without .rows column
#' @param x data.frame
#' @return data.frame without the .rows column
#' @noRd
get_group_data <- function(df){
  grps <- attr(df, "groups")
  grps[, -ncol(grps)]
}
