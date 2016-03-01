#' Map signals over intervals
#' 
#' @param interval_tbl tbl of intervals 
#' @param singal_tbl tbl of signals 
#' @param signal_col bare column name in \code{signal_tbl} for \code{operation}
#' @param operation operation to perform on intersected intervals. One of mean, median,
#'        sum, min, max, absmin, absmax 
#'        
#' @return \code{data.frame}
#' 
#' @examples
#' 
#' signal <- read_bedgraph('inst/extdata/test.bg.gz')
#' intervals <- read_bed('inst/extdata/3fields.bed.gz')
#'  
#' @export
bed_map <- function(interval_tbl, signal_tbl, signal_col,
                    operation = op_choices) {
 
  operation <- match.arg(operation, op_choices) 
  
  intersect_res <- bed_intersect_(a, b, full = TRUE)
  
  map_result <- intersection_result %>%
    group_by(chrom, .start_a, .start_b) %>%
    summarize(result = operation(signal_col)) %>%
    select(chrom, col_spec, result)
  
  map_result 
}

op_choices <- c('sum', 'mean', 'median', 'max', 'min',
                'absmax', 'absmin', 'collapse', 'distinct',
                'count_distinct', 'first', 'last')
