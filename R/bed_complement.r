#' Identify intervals in a genome not covered by a query.
#' 
#' @param x tbl of intervals
#' @param genome chrom sizes
#' 
#' @return \code{data_frame}
#' 
#' @examples 
#' genome <- tibble::frame_data(
#'    ~chrom,  ~size,
#'    "chr1", 500,
#'    "chr2", 600,
#'    "chr3", 800
#' ) 
#' 
#' x <- tibble::frame_data(
#'    ~chrom, ~start, ~end,
#'    "chr1", 100,    300,
#'    "chr1", 200,    400,
#'    "chr2", 1,      100,
#'    "chr2", 200,    400,
#'    "chr3", 500,    600
#' )
#' 
#' # intervals not covered by x
#' bed_complement(x, genome)
#' 
#' @export
bed_complement <- function(x, genome) {

  if ( ! is_merged(x) ) {
    res <- bed_merge(x)
  } 

  # tbl is sorted at this point
  lags <- x %>% group_by(chrom) %>% mutate(.prev_end = lag(end))
  
  first <- lags %>% filter(is.na(.prev_end) & start > 1) %>%
    mutate(.start = 1, .end = start) %>%
    mutate(start = .start, end = .end) %>%
    select(-.prev_end, -.start, -.end)
  
  internal <- lags %>%
    filter(!is.na(.prev_end)) %>%
    mutate(.start = .prev_end, .end = start) %>%
    mutate(start = .start, end = .end) %>%
    select(-.prev_end, -.start, -.end)
  
  final <- lags %>%
    summarize(max.end = max(end)) %>%
    left_join(genome, by = 'chrom') %>%
    filter(size != max.end) %>%
    mutate(start = max.end, end = size) %>%
    select(-size, -max.end)
  
  res <- bind_rows(list(first, internal, final))
  
  res <- bed_sort(res)
 
  res 
}
