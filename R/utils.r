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
#' \code{format_bed} returns a tbl whose columns are ordered by another tbl. 
#' The \code{x} tbl columns are reordered based on the \code{y} columns ordering.
#' If there are \code{x} columns that do not exist in \code{y} they are moved to the last column. 
#'
#'  
#' @param x tbl of intervals
#' @param y tbl of intervals
#'  
#' @examples
#' x <- tibble::tribble(
#'   ~end,  ~chrom,   ~start, ~value,
#'   75,  "chr1",    125,    10
#'  )
#' 
#' y <- tibble::tribble(
#'   ~chrom,   ~start,    ~end,  ~scores,
#'   "chr1",    50,       100,  1.2,
#'   "chr1",    100,       150,  2.4
#'   )
#'   
#' 
#' format_bed(x, y)
#' @noRd
format_bed <- function(x, y) {
  names_x <- names(x)
  names_y <- names(y)
  
  names_x <- names_x[order(match(names_x,names_y))]
  
  x <- select(x, one_of(names_x))
  x
}


#' Identify groups shared between to tbls
#'
#' Identify minimum shared groups between \code{x} and \code{y} tbls. Returns
#' \code{NULL} if there are no shared groups.
#'
#' @param x tbl of intervals
#' @param y tbl of intervals
#'  
#' @return \code{list} of groups or \code{NULL}
#'  
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 150,    400,  100,
#'   "chr2", 230,    430,  200
#' )
#' 
#' y <- tibble::tribble(
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
