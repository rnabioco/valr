#https://github.com/arq5x/bedtools2/blob/master/test/closest/test-closest.sh
context("bed_closest")


test_that("1bp closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  ) 
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal test of 1bp closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("0bp apart closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)   
})

test_that("reciprocal of 0bp apart closer, check for off-by-one errors", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1",	20, 21
  ) 
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)   
})

test_that("check that first left interval at index 0 is not lost", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 9, 10
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that first right interval at index 0 is not lost", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1",	10,	20
  )
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end,
    "chr1", 20, 21
  ) 
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1) 
}
)

test_that("check that strand overlap works", {
  x <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10,	"+")
  
  y <- tibble::frame_data(
    ~chrom,   ~start,    ~end, ~name, ~score, ~strand,
    "chr1", 90, 120, "b", 1,	"-")
  
  res <- bed_closest(x, y, strand = TRUE)
  expect_equal(nrow(res), 0)
}
)
  
