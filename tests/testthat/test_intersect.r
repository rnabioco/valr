context('intersect_cpp')

test_that("simple overlap", {
   bed1_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    150,       250,
    "chr1",    400,       500
  )
  
  bed2_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    175,       200,
    "chr1",    175,       225
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 4)   
})

test_that("multple a's", {
   bed1_df <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200
  )
  
  bed2_df <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    175,       200
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 5)   
})

test_that("multple b's", {
   bed1_df <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200
  )
  
  bed2_df <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200
  )
  
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_equal(nrow(res), 5)   
})

test_that("no overlaps returns empty df", {
  bed1_df <- tibble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  )
  bed2_df <- tibble(
    ~chrom, ~start, ~end,
    "chr1", 300,    400
  )
  res <- intersect_cpp(bed1_df, bed2_df)
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 0)
}) 
  
