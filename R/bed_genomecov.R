#' Calculate coverage across a genome
#'
#' This function is useful for calculating interval coverage across an entire genome.
#'
#' @param x [ivl_df]
#' @param genome [genome_df]
#' @param zero_depth If TRUE, report intervals with zero depth. Zero depth intervals will
#' be reported with respect to groups.
#'
#' @template groups
#'
#' @return
#' [ivl_df] with the an additional column:
#'
#'   - `.depth` depth of interval coverage
#'
#' @family single set operations
#'
#' @seealso
#' \url{https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html}
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~strand,
#'   "chr1", 20, 70, "+",
#'   "chr1", 50, 100, "-",
#'   "chr1", 200, 250, "+",
#'   "chr1", 220, 250, "+"
#' )
#'
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 500,
#'   "chr2", 1000
#' )
#'
#' bed_genomecov(x, genome)
#'
#' bed_genomecov(dplyr::group_by(x, strand), genome)
#'
#' bed_genomecov(dplyr::group_by(x, strand), genome, zero_depth = TRUE)
#'
#' @export
bed_genomecov <- function(x, genome, zero_depth = FALSE) {
  check_required(x)
  check_required(genome)

  x <- check_interval(x)
  genome <- check_genome(genome)

  non_genome_chroms <- setdiff(unique(x$chrom), genome$chrom)
  if (length(non_genome_chroms) > 0) {
    cli::cli_warn(c(
      paste0(
        "The following chromosomes in bed intervals are not",
        " in the genome and will be ignored:"
      ),
      paste0(
        non_genome_chroms,
        sep = "\n"
      )
    ))
    x <- x[!x[["chrom"]] %in% non_genome_chroms, ]
  }

  grp_cols <- group_vars(x)
  x <- bed_sort(x)

  groups <- rlang::syms(unique(c("chrom", grp_cols)))
  x <- group_by(x, !!!groups)

  max_coords <- group_chrom_sizes(x, genome)

  res <- gcoverage_impl(x, max_coords)
  res <- tibble::as_tibble(res)

  # drop non-grouped cols as values no longer match ivls
  res <- select(res, chrom, start, end, one_of(grp_cols), .depth)

  if (!zero_depth) {
    res <- res[res[[".depth"]] > 0, ]
  } else {
    # handle any missing chromosome, zero depth intervals handled on cpp side
    missing_chroms <- setdiff(genome$chrom, group_data(x)$chrom)
    if (length(missing_chroms) > 0) {
      missing_chrom_ivls <- genome[genome[["chrom"]] %in% missing_chroms, ]
      missing_chrom_ivls[["start"]] <- 0L
      missing_chrom_ivls[[".depth"]] <- 0L
      missing_chrom_ivls <- select(missing_chrom_ivls, chrom, start, end = size, .depth)

      if (length(groups) > 1) {
        missing_chrom_ivls <- fill_missing_grouping(missing_chrom_ivls, x)
      }
      res <- bind_rows(res, missing_chrom_ivls)
      res <- bed_sort(res)
    }
  }

  res
}

# return chrom sizes for each group
group_chrom_sizes <- function(x, genome) {
  xx <- left_join(group_data(x), genome, by = "chrom")
  xx[["size"]]
}

# expand df to contain non-chrom groups from grp_df
fill_missing_grouping <- function(df, grp_df) {
  grps <- get_group_data(grp_df)
  grps <- grps[, setdiff(colnames(grps), "chrom")]

  # expand df rows for new groups
  df_grown <- df[rep(seq_len(nrow(df)), nrow(grps)), ]
  # expand groups df for new groups
  grp_grown <- grps[rep(seq_len(nrow(grps)), nrow(df)), ]

  stopifnot(nrow(df_grown) == nrow(grp_grown))
  res <- bind_cols(df_grown, grp_grown)
  res
}
