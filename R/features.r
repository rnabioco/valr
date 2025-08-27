#' Create intron features.
#'
#' Numbers in the `score` column are intron numbers from 5' to 3' independent of
#' strand. I.e., the first introns for `+` and `-` strand genes both have `score`
#' values of `1`.
#'
#' @param x [ivl_df] in BED12 format
#'
#' @family feature functions
#'
#' @examples
#' x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' create_introns(x)
#'
#' @export
create_introns <- function(x) {
  res <- bed12_to_exons(x)
  res <- group_by(res, .data[["name"]])
  res <- mutate(
    res,
    .start = .data[["end"]],
    .end = lead(.data[["start"]]),
    score = ifelse(
      .data[["strand"]] == "+",
      .data[["score"]],
      .data[["score"]] - 1
    ),
    start = .data[[".start"]],
    end = .data[[".end"]]
  )
  res <- select(res, -all_of(c(".start", ".end")))
  res <- ungroup(res)
  res <- na.omit(res)

  # remove zero length intervals
  res <- filter(res, .data[["start"]] < .data[["end"]])

  res
}

#' Create 5' UTR features.
#'
#' @param x [ivl_df] in BED12 format
#'
#' @family feature functions
#'
#' @examples
#' x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' create_utrs5(x)
#'
#' @export
create_utrs5 <- function(x) {
  res <- group_by(x, .data[["name"]])
  res <- mutate(
    res,
    start = ifelse(
      .data[["strand"]] == "+",
      .data[["start"]],
      .data[["cds_end"]]
    ),
    end = ifelse(.data[["strand"]] == "+", .data[["cds_start"]], .data[["end"]])
  )
  res <- ungroup(res)
  res <- select(res, all_of("chrom"):all_of("strand"))

  # remove zero length intervals
  res <- filter(res, .data[["start"]] < .data[["end"]])

  res
}

#' Create 3' UTR features.
#'
#' @param x [ivl_df] in BED12 format
#'
#' @family feature functions
#'
#' @examples
#' x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' create_utrs3(x)
#'
#' @export
create_utrs3 <- function(x) {
  res <- group_by(x, .data[["name"]])
  res <- mutate(
    res,
    start = ifelse(
      .data[["strand"]] == "+",
      .data[["cds_end"]],
      .data[["start"]]
    ),
    end = ifelse(.data[["strand"]] == "+", .data[["end"]], .data[["cds_start"]])
  )
  res <- ungroup(res)
  res <- select(res, all_of("chrom"):all_of("strand"))

  # remove zero length intervals
  res <- filter(res, .data[["start"]] < .data[["end"]])

  res
}

#' Create transcription start site features.
#'
#' @param x [ivl_df] in BED format
#'
#' @family feature functions
#'
#' @examples
#' x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
#'
#' create_tss(x)
#'
#' @export
create_tss <- function(x) {
  res <- group_by(x, .data[["name"]])
  res <- mutate(
    res,
    start = ifelse(
      .data[["strand"]] == "+",
      .data[["start"]],
      .data[["end"]] - 1
    ),
    end = ifelse(.data[["strand"]] == "+", .data[["start"]] + 1, .data[["end"]])
  )
  res <- ungroup(res)
  res <- select(res, all_of("chrom"):all_of("strand"))
  res
}
