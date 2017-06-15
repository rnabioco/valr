context('bed_map')

#test_that("`chrom` grouping throws an error", {
#  x <- tibble::tribble(
#    ~chrom, ~start, ~end,
#    "chr1", 100, 250,
#    "chr2", 250, 500
#  ) %>% group_by(chrom)
#
#  y <- tibble::tribble(
#    ~chrom, ~start, ~end, ~value,
#    "chr1", 100, 250, 10,
#    "chr1", 150, 250, 20,
#    "chr2", 250, 500, 500
#  ) %>% group_by(chrom)
#
#  expect_error(bed_map(x, y))
#})

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

  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~id, ~vals,
    "chr1", 100, 250, 1,  10,
    "chr2", 250, 500, 3,  500,
    "chr2", 250, 500, 2, NA
  )
  res <- bed_map(x, y, vals = sum(value))
  expect_true(all(res == pred, na.rm = TRUE))
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

  res <- bed_map(x, y, vals = values_unique(value))
  expect_equal(res$vals, c("10,20"))
})

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~id,
  "chr1", 100, 200, 1,
  "chr1", 250, 500, 2,
  "chr2", 250, 500, 3
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value,
  "chr1", 100, 150, 10,
  "chr1", 150, 250, 20,
  "chr1", 140, 250, 30,
  "chr1", 150, 200, 40
)

test_that("concat works correctly", {
  res <- bed_map(x, y, vals = concat(value))
  expected <- c('10,30,20,40', NA, NA)
  expect_equal(res$vals, expected)
})

test_that("values works correctly", {
  res <- bed_map(x, y, vals = values(value))
  expected <- c('10,30,20,40', NA, NA)
  expect_equal(res$vals, expected)
})

test_that("first works correctly", {
  res <- bed_map(x, y, first = first(value))
  expected <- c(10, NA, NA)
  expect_equal(res$first, expected)
})

test_that("last works correctly", {
  res <- bed_map(x, y, last = last(value))
  expected <- c(40, NA, NA)
  expect_equal(res$last, expected)
})

test_that("book-ended intervals are not reported", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 100, 150, 10,
    "chr1", 200, 250, 20
  )

  expected <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    'chr1',  100,    200,   10
  )
  res <- bed_map(x, y, value = sum(value))
  expect_equal(res, expected)
})

test_that("ensure that mapping is calculated with respect to input tbls issue#108",{

  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    'chr1', 100,    200,  'B',
    'chr1', 200,    400,  'A',
    'chr1', 500,    600,  'C',
    'chr2', 125,    175,  'C',
    'chr2', 150,    200,  'A',
    'chr3', 100,    300,  'A'
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~value,
    'chr1', 100,    199,  'A',  10,
    'chr1', 200,    400,  'B', 20,
    'chr1', 500,    600,  'A', 30,
    'chr2', 125,    175,  'C', 40,
    'chr2', 350,    500,  'A', 50,
    'chr3', 500,    600,  'A', 100
  )

  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~total,
    'chr1', 100,    200,  'B', NA,
    'chr1', 200,    400,  'A', NA,
    'chr1', 500,    600,  'C', NA,
    'chr2', 125,    175,  'C', 40,
    'chr2', 150,    200,  'A', NA,
    'chr3', 100,    300,  'A', NA
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group)

  res <- bed_map(x, y, total = sum(value))
  expect_true(all(pred == res, na.rm = T))
})

