context('bed_intersect')
 
test_that("simple overlap", {
   x <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    150,       250,
    "chr1",    400,       500
  )
  
  y <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    175,       200,
    "chr1",    175,       225
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 4)   
})

test_that("multple a's", {
   x <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200,
    "chr1",    100,       200
  )
  
  y <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    175,       200
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)   
})

test_that("multple b's", {
   x <- tibble(
    ~chrom,    ~start,    ~end,
    "chr1",    100,       200
  )
  
  y <- tibble(
    ~chrom,   ~start,    ~end,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200,
    "chr1",    175,       200
  )
  
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)   
})

test_that("no overlaps returns empty df", {
  x <- tibble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  )
  y <- tibble(
    ~chrom, ~start, ~end,
    "chr1", 300,    400
  )
  res <- bed_intersect(x, y)
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 0)
}) 
