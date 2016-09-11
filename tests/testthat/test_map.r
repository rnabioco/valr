context('bed_map')

test_that("`chrom` grouping throws an error", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250,
    "chr2", 250, 500
  ) %>% group_by(chrom)
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr2", 250, 500, 500
  ) %>% group_by(chrom)
  
  expect_error(bed_map(x, y))
})

test_that("x/y groupings are respected", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~id,
    "chr1", 100, 250, 1,
    "chr2", 250, 500, 2,
    "chr2", 250, 500, 3
  ) %>% group_by(id)
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value, ~id,
    "chr1", 100, 250, 10,  1,
    "chr1", 150, 250, 20,  2,
    "chr2", 250, 500, 500, 3
  ) %>% group_by(id)
  
  res <- bed_map(x, y, vals = sum(value.y))
  expect_equal(res$vals, c(10,20,500,500))
})

test_that("concat works correctly", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250,
    "chr2", 250, 500,
    "chr2", 250, 500
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr2", 250, 500, 500
  )
  
  res <- bed_map(x, y, vals = concat(value.y))
  expect_equal(res$vals, c("10,20", "500,500"))
})

test_that("values works correctly", {
   x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250
   )
  
   y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 100, 250, 20,
    "chr1", 150, 250, 10,
    "chr1", 150, 250, 20
   )

  res <- bed_map(x, y, vals = values(value.y))
  expect_equal(res$vals, c("10,20,10,20"))
})

test_that("values_unique works correctly", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20
  )
  
  res <- bed_map(x, y, vals = values_unique(value.y))
  expect_equal(res$vals, c("10,20"))
})

test_that("first works correctly", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr1", 100, 250, 30,
    "chr1", 150, 250, 40
  )
  
  res <- bed_map(x, y, vals = first(value.y))
  expect_equal(res$vals, 10)
})

test_that("last works correctly", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 250
  )
  
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr1", 100, 250, 30,
    "chr1", 150, 250, 40
  )
  
  res <- bed_map(x, y, vals = last(value.y))
  expect_equal(res$vals, 40)
})
