context('bed_sort')

test_that('intervals can be sorted by size', {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr1', 200,    400,
    'chr2', 100,    200,
    'chr1', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_size = TRUE))
  expect_equal(res$start, c(100,200,400,1000))
})

test_that('intervals can be sorted by size and chrom', {
    x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr1', 200,    400,
    'chr2', 100,    200,
    'chr1', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_size = TRUE, by_chrom = TRUE))
  expect_equal(res$start, c(100,200,400,1000))
})

test_that('intervals can be reverse sorted by size', {
    x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr1', 200,    400,
    'chr2', 100,    200,
    'chr1', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_size = TRUE, reverse = TRUE))
  expect_equal(res$start, c(1000,400,200,100))
})

test_that('intervals can be reverse sorted by size and chrom', {
    x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr1', 200,    400,
    'chr2', 100,    200,
    'chr1', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_size = TRUE, by_chrom = TRUE, reverse = TRUE))
  expect_equal(res$start, c(1000,400,200,100))
})

test_that('intervals can be sorted by chrom', {
    x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr2', 100,    200,
    'chr1', 200,    400,
    'chr3', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_chrom = TRUE))
  expect_equal(res$start, c(200,1000,100,400))
})

test_that('intervals can be reverse sorted by start and chrom', {
   x <- tibble::tribble(
    ~chrom, ~start, ~end,
    'chr1', 1000,   2000,
    'chr2', 100,    200,
    'chr1', 200,    400,
    'chr3', 400,    700
    )
  expect_warning(res <- bed_sort(x, by_chrom = TRUE, reverse = TRUE))
  expect_equal(res$start, c(1000,200,100,400)) 
})
