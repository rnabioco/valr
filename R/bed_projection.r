#' Projection Test for query interval overlap.
#' 
#' @param x tbl of intervals
#' @param y tbl of intervals
#' @param genome chrom sizes
#' @param by_chrom compute test per chromosome  rather than across whole genome
#' 
#' @template groups
#' 
#' @family interval-stats
#' @return \code{data_frame} with the following columns:
#'   \itemize{ 
#'     \item{\code{chrom}} {the name of chromosome tested if \code{by_chrom} is \code{TRUE}, otherwise set to \code{whole_genome}}
#'     \item{\code{.p_value}} {p-value from a binomial test, note that p-values > 0.5 will be reported as 1 - p-value and .lower_tail will be set to FALSE}
#'     \item{\code{.effect_ratio}} {ratio of observed overlap probabilities to expected}
#'     \item{\code{.lower_tail}} {TRUE denotes thate observed probabilities are in the lower tail of the distribution (less overlap than expected), FALSE
#'     denotes that the observed probabililty are in the upper tail of the distribution (more overlap than expected)}
#'     }
#'   
#' @seealso
#'   \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'   
#' @examples 
#' genome <- tibble::tribble(
#'  ~chrom, ~size,
#'  "chr1", 1e6,
#'  "chr2", 2e6,
#'  "chr3", 4e6
#' )
#' 
#' x <- bed_random(genome)
#' bed_shuffle(x, genome)
#' 
#' @export
bed_projection <- function(x, y, genome, by_chrom = FALSE) {

  #find midpoints
  x <- mutate(x, 
              .midpoint = round((end + start) / 2),
              start = .midpoint,
              end = .midpoint + 1)
  x <- select(x, -.midpoint)

  # count overlaps per chromosome
  obs_counts <- bed_intersect(x, y)
  
  # count overlaps
  obs_counts <- group_by(obs_counts, chrom)
  obs_counts <- summarize(obs_counts, .obs_counts = n())
  
  # total x intervals tested
  obs_counts <- mutate(obs_counts, .total_trials = nrow(x))
  
  # calcuate probabilty of overlap by chance
  y <- mutate(y, 
              .length = end - start)
  y <- group_by(y, chrom)
  y <- summarize(y, .reference_coverage = sum(.length))
  y <- inner_join(y, genome, by = "chrom")
  null_dist <- mutate(y, .exp_prob = .reference_coverage / size)
 
  res <- inner_join(obs_counts, null_dist, by = "chrom")
  
  # binomial test and obs/exp
  if (by_chrom){
    res <- group_by(res, chrom)
    res <- summarize(res,
                     .p_value = pbinom(q = .obs_counts, 
                                      size = .total_trials,
                                      prob = .exp_prob),
                     .effect_ratio = (.obs_counts / .total_trials) / .exp_prob
                     )
  } else {
    res <- ungroup(res)
    res <- summarize(res,
                     .obs_counts = sum(.obs_counts),
                     .total_trials = sum(.total_trials),
                     .exp_prob = sum(.reference_coverage) / sum(size))
    
    res <- summarize(res,
                     chrom = "whole_genome",
                     .p_value = pbinom(q = .obs_counts, 
                                      size = .total_trials,
                                      prob = .exp_prob),
                     .effect_ratio = (.obs_counts / .total_trials) / .exp_prob)
  }
  res <- mutate(res,
                .lower_tail = if_else(.p_value < .5,
                                      "TRUE",
                                      "FALSE"),
                .p_value = if_else(.p_value < .5,
                                      .p_value,
                                      1 - .p_value))
  res
}  













