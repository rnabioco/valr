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

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~id,
  "chr1", 100, 200, 1,
  "chr1", 250, 500, 2,
  "chr2", 250, 500, 3
) %>% group_by(id)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value,
  "chr1", 100, 150, 10,
  "chr1", 150, 250, 20,
  "chr1", 140, 250, 30,
  "chr1", 150, 200, 40
)

test_that("concat works correctly", {
  res <- bed_map(x, y, vals = concat(value.y))
  expect_equal(res$vals, c("10,30,20,40", "30,20"))
})

test_that("values works correctly", {
  res <- bed_map(x, y, vals = values(value.y))
  expect_equal(res$vals, c("10,30,20,40", "30,20"))
})

test_that("first works correctly", {
  res <- bed_map(x, y, first = first(value.y))
  expect_true(all(res$first == c(10, 30)))
})

test_that("last works correctly", {
  res <- bed_map(x, y, last = last(value.y))
  expect_true(all(res$last == c(40, 20))) 
})




