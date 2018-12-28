#' Identify intervals in a genome not covered by a query.
#'
#' @param x [tbl_interval()]
#' @param genome [tbl_genome()]
#'
#' @family single set operations
#'
#' @return [tbl_interval()]
#'
#' @examples
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 0,      10,
#'   'chr1', 75,     100
#' )
#'
#' genome <- trbl_genome(
#'   ~chrom, ~size,
#'   'chr1', 200
#' )
#'
#' bed_glyph(bed_complement(x, genome))
#'
#' genome <- trbl_genome(
#'    ~chrom,  ~size,
#'    'chr1',  500,
#'    'chr2',  600,
#'    'chr3',  800
#' )
#'
#' x <- trbl_interval(
#'    ~chrom, ~start, ~end,
#'    'chr1', 100,    300,
#'    'chr1', 200,    400,
#'    'chr2', 0,      100,
#'    'chr2', 200,    400,
#'    'chr3', 500,    600
#' )
#'
#' # intervals not covered by x
#' bed_complement(x, genome)
#'
#' @export
bed_complement <- function(x, genome) {
  if (!is.tbl_interval(x)) x <- as.tbl_interval(x)
  if (!is.tbl_genome(genome)) genome <- as.tbl_genome(genome)

  res <- bed_merge(x)

  # non-overlapping chroms
  chroms_no_overlaps <- anti_join(genome, res, by = "chrom")
  chroms_no_overlaps <- mutate(chroms_no_overlaps, start = 0)
  chroms_no_overlaps <- select(chroms_no_overlaps, chrom, start, end = size)

  # remove rows from x that are not in genome
  res <- semi_join(res, genome, by = "chrom")

  res <- group_by(res, chrom)

  if (utils::packageVersion("dplyr") < "0.7.99.9000"){
    res <- update_groups(res)
  }

  res <- complement_impl(res, genome)

  res <- bind_rows(res, chroms_no_overlaps)
  res <- bed_sort(res)

  res
}
