context('bed_random')

genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 1e6
) 

test_that("random takes genome on input", {
  res <- bed_random(genome, n = 1e6)
  expect_equal(nrow(res), 1e6)   
})
