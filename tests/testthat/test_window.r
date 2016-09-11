context('bed_window')

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 5000,
  "chr2", 10000
)

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 500,    1000, '.',   '.',     '+',
  "chr1", 1000,   1500, '.',   '.',     '-',
  "chr2", 1000,   1200, '.',   '.',     '-'
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 400,    450, '.',   '.',     '+',
  "chr1", 1000,   1200, '.',   '.',     '-',
  "chr1", 1100,    1500, '.',   '.',     '+',
  "chr2", 1300,   1500, '.',   '.',     '-'
)

test_that("both arg works", {
  out <- bed_window(x, y, genome, both = 110)
  expect_equal(nrow(out), 6)
})

test_that("left arg works", {
  out <- bed_window(x, y, genome, left = 110)
  expect_equal(nrow(out), 4)
})

test_that("right arg works", {
  out <- bed_window(x, y, genome, right = 110)
  expect_equal(nrow(out), 5)
})

test_that("strand arg works", {
  out <- bed_window(x, y, genome, right = 110, strand = TRUE)
  expect_equal(nrow(out), 2)
})

test_that("strand_opp arg works", {
  dist <- 100
  out <- bed_window(x, y, genome, right = 110, 
                    strand_opp = TRUE, strand = TRUE)

  expect_equal(nrow(out), 2)

})

# test fraction
