#' Projection test for query interval overlap.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome chrom sizes
#' @param by_chrom compute test per chromosome
#' 
#' @template stats
#' 
#' @family interval-stats
#' 
#' @return \code{data_frame} with the following columns:
#'   \itemize{ 
#'     \item{\code{chrom}} {the name of chromosome tested if \code{by_chrom = TRUE},
#'      otherwise has a value of \code{'whole_genome'}}
#'     \item{\code{p.value}} {p-value from a binomial test. p-values > 0.5
#'      will be reported as 1 - p-value and \code{lower_tail} will be \code{FALSE}}
#'     \item{\code{obs_exp_ratio}} {ratio of observed to expected overlap frequency}
#'     \item{\code{lower_tail}} a boolean column. {\code{TRUE} indicates the observed overlaps
#'      is in the lower tail of the distribution (e.g., less overlap than expected); \code{FALSE}
#'      indicates that the observed overlaps are in the upper tail of the distribution
#'      (e.g., more overlap than expected)}
#'     }
#'   
#' @seealso
#'   \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'   
#' @examples 
#' genome <- tibble::tribble(
#'  ~chrom, ~size,
#'  "chr1", 1e4,
#'  "chr2", 2e4,
#'  "chr3", 4e4
#' )
#' 
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 100,    200,
#'  "chr1", 250,    400,
#'  "chr1", 500,    600,
#'  "chr1", 1000,   2000,
#'  "chr2", 100,    200
#' )
#'
#' y <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  "chr1", 150,    175,
#'  "chr1", 525,    575,
#'  "chr1", 1100,   1200,
#'  "chr1", 1400,   1600,
#'  "chr2", 200,    1500
#' ) 
#'
#' bed_projection(x, y, genome)
#' bed_projection(x, y, genome, by_chrom = TRUE)
#' 
#' @export
bed_projection <- function(x, y, genome, by_chrom = FALSE) {

  #find midpoints
  x <- mutate(x, .midpoint = round((end + start) / 2),
              start = .midpoint,
              end = .midpoint + 1)
  x <- select(x, -.midpoint)
  
  # flatten y intervals
  y <- bed_merge(y)
  
  # count overlaps per chromosome, 
  obs_counts <- bed_intersect(x, y)

  # count overlaps
  obs_counts <- group_by(obs_counts, chrom)
  obs_counts <- summarize(obs_counts, .obs_counts = n())
          
  #total x intervals tested
  total_counts <- group_by(x, chrom)
  total_counts <- summarize(total_counts, .total_trials = n())
  obs_counts <- full_join(obs_counts, total_counts, by = "chrom")
  obs_counts <- mutate(obs_counts, .obs_counts = if_else(is.na(.obs_counts), 
                                                         as.integer(0), .obs_counts))
  
  # calculate probabilty of overlap by chance
  y <- mutate(y, .length = end - start)
  y <- group_by(y, chrom)
  y <- summarize(y, .reference_coverage = sum(.length))
  
  # add in any missing chromosomes
  y <- full_join(y, genome, by = "chrom")
  
  null_dist <- mutate(y, .exp_prob = .reference_coverage / size)
 
  res <- inner_join(obs_counts, null_dist, by = "chrom")
  
  # binomial test and obs/exp
  if (by_chrom){
    res <- group_by(res, chrom)
    res <- summarize(res,
                     p.value = stats::pbinom(q = .obs_counts, 
                                      size = .total_trials,
                                      prob = .exp_prob),
                     obs_exp_ratio = (.obs_counts / .total_trials) / .exp_prob
                     )
  } else {
    res <- ungroup(res)
    res <- summarize(res,
                     .obs_counts = sum(.obs_counts),
                     .total_trials = sum(.total_trials),
                     .exp_prob = sum(.reference_coverage) / sum(size))
    
    res <- summarize(res,
                     chrom = "whole_genome",
                     p.value = stats::pbinom(q = .obs_counts, 
                                      size = .total_trials,
                                      prob = .exp_prob),
                     obs_exp_ratio = (.obs_counts / .total_trials) / .exp_prob)
  }
  
  res <- mutate(res,
                lower_tail = if_else(p.value < .5,
                                     "TRUE",
                                     "FALSE"),
                p.value = if_else(p.value < .5,
                                  p.value,
                                  1 - p.value))
  res
}  
