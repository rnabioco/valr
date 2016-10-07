#' provides working directory for valr example files
#' 
#' @export
valr_example <- function(path) {
  # https://twitter.com/JennyBryan/status/780150538654527488
  system.file("extdata", path, package = 'valr', mustWork = TRUE)
}

#' reformat bed tbl to match another tbl
#' 
#' \code{format_bed} returns a tbl whose columns are ordered by another tbl.
#'Columns not found in y tbl are dropped. Missing y columns are added and populated
#' with a dummy entry "."
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

}

