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

test_that("merge with interval at start", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,         50,
    "chr1",    100,       200,
    "chr1",    150,       250
  )
  
  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 2)   
})

test_that("merge with two chroms", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,         50,
    "chr1",    25,        75,
    "chr2",    100,       200,
    "chr2",    150,       250
  )
  
  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 2)   
})

test_that("book-ended intervals are no merged", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,         50,
    "chr1",    50,        100
  )
  
  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 2)   
})

test_that("max_dist is enforced", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,         50,
    "chr1",    50,        100
  )
  
  res <- bed_merge(bed_df, max_dist = 1)
  expect_equal(nrow(res), 1)   
})

test_that("max_dist is a positive value", {
  bed_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    1,         50,
    "chr1",    50,        100
  )
  
  expect_error(bed_merge(bed_df, max_dist = -1))
})
