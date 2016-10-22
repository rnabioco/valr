#' provides working directory for valr example files
#' @param path path to file
#' @export
valr_example <- function(path) {
  # https://twitter.com/JennyBryan/status/780150538654527488
  system.file("extdata", path, package = 'valr', mustWork = TRUE)
}

#' reformat bed tbl to match another tbl
#' 
#' \code{format_bed} returns a tbl whose columns are ordered by another tbl.
#' Columns not found in Y tbl are dropped from X. Y columns not found in X 
#' added to X and populated with a dummy entry "."
#' @param x tbl of intervals
#' @param y tbl of intervals 
#' @export
format_bed <- function(x, y) {
  names_x <- names(x)
  names_y <- names(y)
  
  if (any(!names_y %in% names_x)){
    cols_to_add <- setdiff(names_y, names_x)
    n <- ncol(x)
    for (i in seq_along(cols_to_add)){
      x[n + i] <- "."
      colnames(x)[n + i] <- cols_to_add[i]
    }
  }
  x <- select(x, one_of(names_y))
  x
}


#' Compare two tbl_dfs and find the minimum shared groups 
#' 
#' \code{shared_groups} returns the minimum shared groups between
#' \code{x} and \code{y} tbls. Return value is \code{NULL} if there are
#' no shared groups. 
#' 
#' @param x tbl 
#' @param y tbl 
#' @return \code{list}
#' @export

shared_groups <- function(x, y) {
  groups_x <- groups(x)
  groups_y <- groups(y)
  
  groups_xy <- intersect(groups_x, groups_y)
  if (length(groups_xy) == 0){
    groups_xy <- NULL
  }
  groups_xy
}

