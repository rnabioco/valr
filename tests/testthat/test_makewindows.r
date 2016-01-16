context("bed_makewindows")

genome <- dplyr::tibble(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 400
)

bed_df <- dplyr::tibble(
  ~chrom, ~start, ~end, ~name,
  "chr1", 100, 200, 'A',
  "chr2", 300, 350, 'B'
) 

res <- bed_makewindows(bed_df, genome, num_windows = 10)
expect_that(nrow(res), 10)