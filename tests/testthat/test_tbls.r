context("tbls")

test_that("tbl_interval classes are set", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1',  1,     50
  )
  x <- tbl_interval(x)
  expect_true(is.tbl_interval(x))
})

test_that("tbl_genome classes are set", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1',  1e4
  )
  genome <- tbl_genome(genome)
  expect_true(is.tbl_genome(genome))
})

test_that("invalid tbl_interval names throw error", {
  x <- tibble::tribble(
    ~pork, ~pie, ~hat,
    'chr1',  1,  50
  )
  expect_error(tbl_interval(x), 'expected 3 required names, missing: chrom, start, end')

  # missing 1 only
  x <- tibble::tribble(
    ~chrom, ~start, ~oops,
    'chr1',  1,  50
  )
  expect_error(tbl_interval(x), 'expected 3 required names, missing: end')
})

test_that("invalid tbl_genome names throw error", {
  genome <- tibble::tribble(
    ~foo,   ~bar,
    'chr1', 1e4
  )
  expect_error(tbl_genome(genome), 'expected 2 required names, missing: chrom, size')
})

test_that("duplicate tbl_genome refs throw error", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    'chr1', 1e4,
    'chr1', 1e4
  )
  expect_error(tbl_genome(genome), 'duplicate chroms in genome: chr1')
})

test_that("trbl_interval accepts tribble format", {
  x <- trbl_interval(
    ~chrom, ~start, ~end,
    'chr1', 1,      100
  )
  expect_true(is.tbl_interval(x))
  expect_error(trbl_interval(1,2,3))
})

test_that("trbl_genome accepts tribble format", {
   x <- trbl_genome(
    ~chrom, ~size,
    'chr1', 1e6
  )
  expect_true(is.tbl_genome(x))
  expect_error(trbl_genome(1,2,3))
})
