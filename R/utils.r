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
#' Columns not found in \code{y} tbl are dropped from \code{x}. \code{y} columns
#'  not found in \code{x} are added to \code{x} and populated with a dummy entry \code{"."}
#'  
#' @param x tbl of intervals
#' @param y tbl of intervals
#'  
#' @examples
#' x <- tibble::frame_data(
#'   ~end,  ~chrom,   ~start, ~value,
#'   75,  "chr1",    125,    10
#'  )
#' 
#' y <- tibble::frame_data(
#'   ~chrom,   ~start,    ~end,  ~scores,
#'   "chr1",    50,       100,  1.2,
#'   "chr1",    100,       150,  2.4
#'   )
#'   
#' 
#' format_bed(x, y)
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
#' 
#' @examples
#'x <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~value,
#'   "chr1", 150,    400,  100,
#'   "chr2", 230,    430,  200,
#'   )
#' 
#' y <- tibble::frame_data(
#'   ~chrom,   ~start,    ~end,  ~value,
#'   "chr1",    50,       100,  1,
#'   "chr1",    100,       150,  2
#'   )
#'   
#' x <- dplyr::group_by(x, chrom, value)
#' y <- dplyr::group_by(y, chrom, value)
#'shared_groups(x, y)
#'  
#' y <- dplyr::group_by(y, chrom)
#' shared_groups(x, y)
#' 
#' y <-dplyr::ungroup(y)
#' shared_groups(x, y)
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

