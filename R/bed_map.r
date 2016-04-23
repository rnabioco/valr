#' Map signals over intervals
#' 
#' @param bed_tbl tbl of intervals 
#' @param signal_tbl tbl of signals 
#' @param signal_col bare column name in \code{signal_tbl} for \code{operation}
#' @param operation operation to perform on intersected intervals. One of mean, median,
#'        sum, min, max, absmin, absmax 
#'        
#' @return \code{data_frame}
#' 
#' @examples
#' bed_tbl <- tibble::frame_data(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100, 250,
#'  "chr2", 250, 500)
#'  
#' signal_tbl <- tibble::frame_data(
#'  ~chrom, ~start, ~end, ~value,
#'  "chr1", 100, 250, 10,
#'  "chr1", 150, 250, 20,
#'  "chr2", 250, 500, 500)
#' 
#' bed_map(bed_tbl, signal_tbl, value, 'sum')
#' 
#' @export
bed_map <- function(bed_tbl, signal_tbl, signal_col, operation) {
 
  operation <- match.arg(operation, op_choices) 
  
  isect_res <- bed_intersect(bed_tbl, signal_tbl)
  
  map_result <- isect_res %>%
    group_by(chrom, start.x, end.x) %>%
    summarize(.result = sum(value)) %>%
    ungroup()
 
  colnames(map_result) <- str_replace(colnames(map_result), '.x$', '') 
  
  map_result 
}

op_choices <- c('sum', 'mean', 'median', 'max', 'min',
                'absmax', 'absmin', 'collapse', 'distinct',
                'count_distinct', 'first', 'last')
