#' Projection test for query interval overlap.
#'
#' @param x [ivl_df]
#' @param y [ivl_df]
#' @param genome [genome_df]
#' @param by_chrom compute test per chromosome
#'
#' @template stats
#'
#' @family interval statistics
#'
#' @return
#' [ivl_df] with the following columns:
#'
#'   - `chrom` the name of chromosome tested if `by_chrom = TRUE`,
#'      otherwise has a value of `whole_genome`
#'
#'   - `p.value` p-value from a binomial test. p-values > 0.5
#'      are converted to `1 - p-value` and `lower_tail` is `FALSE`
#'
#'   - `obs_exp_ratio` ratio of observed to expected overlap frequency
#'
#'   - `lower_tail` `TRUE` indicates the observed overlaps are in the lower tail
#'     of the distribution (e.g., less overlap than expected). `FALSE` indicates
#'     that the observed overlaps are in the upper tail of the distribution (e.g.,
#'     more overlap than expected)
#'
#' @seealso
#'   \url{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}
#'
#' @examples
#' genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
#'
#' x <- bed_random(genome, seed = 1010486)
#' y <- bed_random(genome, seed = 9203911)
#'
#' bed_projection(x, y, genome)
#'
#' bed_projection(x, y, genome, by_chrom = TRUE)
#'
#' @export
bed_projection <- function(x, y, genome, by_chrom = FALSE) {
  check_required(x)
  check_required(y)
  check_required(genome)

  x <- check_interval(x)
  y <- check_interval(y)
  genome <- check_genome(genome)

  # find midpoints
  x <- mutate(
    x,
    .midpoint = round((.data[["end"]] + .data[["start"]]) / 2),
    start = .data[[".midpoint"]],
    end = .data[[".midpoint"]] + 1
  )
  x <- select(x, -all_of(".midpoint"))

  # flatten y intervals
  y <- bed_merge(y)

  # count overlaps per chromosome,
  # TODO: revisit when min_overlap default changes to 1L
  obs_counts <- bed_intersect(x, y, min_overlap = 0L)

  # count overlaps
  obs_counts <- group_by(obs_counts, .data[["chrom"]])
  obs_counts <- summarize(obs_counts, .obs_counts = n())

  # total x intervals tested
  total_counts <- group_by(x, .data[["chrom"]])
  total_counts <- summarize(total_counts, .total_trials = n())
  obs_counts <- full_join(obs_counts, total_counts, by = "chrom")
  obs_counts <- mutate(
    obs_counts,
    .obs_counts = if_else(
      is.na(.data[[".obs_counts"]]),
      as.integer(0),
      .data[[".obs_counts"]]
    )
  )

  # calculate probabilty of overlap by chance
  y <- mutate(y, .length = .data[["end"]] - .data[["start"]])
  y <- group_by(y, .data[["chrom"]])
  y <- summarize(y, .reference_coverage = sum(.data[[".length"]]))

  # add in any missing chromosomes
  y <- full_join(y, genome, by = "chrom")

  null_dist <- mutate(
    y,
    .exp_prob = .data[[".reference_coverage"]] / .data[["size"]]
  )

  res <- inner_join(obs_counts, null_dist, by = "chrom")

  # binomial test and obs/exp
  if (by_chrom) {
    res <- group_by(res, .data[["chrom"]])
    res <- summarize(
      res,
      p.value = stats::pbinom(
        q = .data[[".obs_counts"]],
        size = .data[[".total_trials"]],
        prob = .data[[".exp_prob"]]
      ),
      obs_exp_ratio = (.data[[".obs_counts"]] / .data[[".total_trials"]]) /
        .data[[".exp_prob"]]
    )
  } else {
    res <- ungroup(res)
    res <- summarize(
      res,
      .obs_counts = sum(.data[[".obs_counts"]]),
      .total_trials = sum(.data[[".total_trials"]]),
      .exp_prob = sum(.data[[".reference_coverage"]]) /
        sum(as.numeric(.data[["size"]]))
    )

    res <- summarize(
      res,
      chrom = "whole_genome",
      p.value = stats::pbinom(
        q = .data[[".obs_counts"]],
        size = .data[[".total_trials"]],
        prob = .data[[".exp_prob"]]
      ),
      obs_exp_ratio = (.data[[".obs_counts"]] / .data[[".total_trials"]]) /
        .data[[".exp_prob"]]
    )
  }

  res <- mutate(
    res,
    lower_tail = if_else(.data[["p.value"]] < .5, "TRUE", "FALSE"),
    p.value = if_else(
      .data[["p.value"]] < .5,
      .data[["p.value"]],
      1 - .data[["p.value"]]
    )
  )
  res
}
