context('bed_merge')

test_that("merge on 1 chrom", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    150,       250,
    "chr1",    200,       350
  )
  
  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 1)   
})

