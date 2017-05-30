#' Create intron features.
#'
#' Numbers in the `score` column are intron numbers from 5' to 3' independent of
#' strand. I.e., the first introns for `+` and `-` strand genes both have `score`
#' values of `1`.
#'
#' @param x [tbl_interval] in BED12 format
#'
#' @examples
#' x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
#'
#' create_introns(x)
#'
#' @export
create_introns <- function(x) {

  res <- bed12_to_exons(x)
  res <- group_by(res, name)
  res <- mutate(res,
                .start = end, .end = lead(start),
                score = ifelse(score == 1, 1, score - 1),
                start = .start, end = .end)
  res <- select(res, -.start, -.end)
  res <- ungroup(res)
  res <- na.omit(res)

  # remove zero length intervals
  res <- filter(res, start < end)

  res
}

#' Create 5' UTR features.
#'
#' @param x [tbl_interval] in BED12 format
#'
#' @examples
#' x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
#'
#' create_utrs5(x)
#'
#' @export
create_utrs5 <- function(x) {
  res <- group_by(x, name)
  res <- mutate(res,
                start = ifelse(strand == '+', start, cds_end),
                end   = ifelse(strand == '+', cds_start, end))
  res <- ungroup(res)
  res <- select(res, chrom:strand)

  # remove zero length intervals
  res <- filter(res, start < end)

  res
}

#' Create 3' UTR features.
#'
#' @param x [tbl_interval] in BED12 format
#'
#' @examples
#' x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
#'
#' create_utrs3(x)
#'
#' @export
create_utrs3 <- function(x) {
  res <- group_by(x, name)
  res <- mutate(res,
                start = ifelse(strand == '+', cds_end, start),
                end   = ifelse(strand == '+', end, cds_start))
  res <- ungroup(res)
  res <- select(res, chrom:strand)

  # remove zero length intervals
  res <- filter(res, start < end)

  res
}
