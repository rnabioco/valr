context("bed_slop")

genome <- tibble::tribble(
  ~ chrom, ~ size,
  "chr1", 5000
)

x <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 500, 1000, ".", ".", "+",
  "chr1", 1000, 1500, ".", ".", "-"
)

test_that("left arg works", {
  dist <- 100
  out <- x %>%
    bed_slop(genome, left = dist)
  expect_true(all(x$start - out$start == dist))
})

test_that("right arg works", {
  dist <- 100
  out <- x %>%
    bed_slop(genome, right = dist)
  expect_true(all(out$end - x$end == dist))
})

test_that("both arg works", {
  dist <- 100
  out <- x %>%
    bed_slop(genome, both = dist)
  expect_true(all(x$start - out$start == dist))
  expect_true(all(out$end - x$end == dist))
})

test_that("both with fraction works", {
  res <- bed_slop(x, genome, both = 0.5, fraction = TRUE)
  expect_equal(res$start, c(250, 750))
  expect_equal(res$end, c(1250, 1750))
})

test_that("left / right with fraction works", {
  res <- bed_slop(x, genome, left = 0.5, fraction = TRUE)
  expect_equal(res$start, c(250, 750))
  expect_equal(res$end, c(1000, 1500))
})

test_that("left, fraction, strand works", {
  res <- bed_slop(x, genome, left = 0.5, fraction = TRUE, strand = TRUE)
  expect_equal(res$start, c(250, 1000))
  expect_equal(res$end, c(1000, 1750))
})

test_that("right, fraction, strand works", {
  res <- bed_slop(x, genome, right = 0.5, fraction = TRUE, strand = TRUE)
  expect_equal(res$start, c(500, 750))
  expect_equal(res$end, c(1250, 1500))
})

test_that("strand with left works", {
  res <- bed_slop(x, genome, left = 100, strand = TRUE)
  expect_equal(res$start, c(400, 1000))
  expect_equal(res$end, c(1000, 1600))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/slop/test-slop.sh
a <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 100L, 200L, "a1", 1L, "+",
  "chr1", 100L, 200L, "a2", 2L, "-"
)

tiny.genome <- tibble::tribble(
  ~ chrom, ~ size,
  "chr1", 1000
)

h19 <- read.table(file = valr_example("hg19.chrom.sizes.gz"), sep = "\t", header = FALSE, stringsAsFactors=FALSE)
colnames(h19) <- c("chrom", "size")
h19 <- tibble::as_tibble(h19)

test_that("test going beyond the start of the chrom", {
  res <- bed_slop(a, tiny.genome, both = 200, trim = TRUE)
  expect_equal(res$start, c(0, 0))
  expect_equal(res$end, c(400, 400))
})

test_that("test going beyond the end of the chrom", {
  res <- bed_slop(a, tiny.genome, left = 0, right = 1000, trim = TRUE)
  expect_equal(res$start, c(100, 100))
  expect_equal(res$end, c(1000, 1000))
})

test_that("test going beyond the start and end of the chrom", {
  res <- bed_slop(a, tiny.genome, both = 2000, trim = TRUE)
  expect_equal(res$start, c(0, 0))
  expect_equal(res$end, c(1000, 1000))
})

test_that("test going beyond the start and end of the chrom with strand", {
  res <- bed_slop(a, tiny.genome, both = 2000, strand = TRUE, trim = TRUE)
  expect_equal(res$start, c(0, 0))
  expect_equal(res$end, c(1000, 1000))
})

test_that("test slop factor being larger than a signed int", {
  res <- bed_slop(a, tiny.genome, both = 3000000000, strand = TRUE, trim = TRUE)
  expect_equal(res$start, c(0, 0))
  expect_equal(res$end, c(1000, 1000))
})

test_that("test that old floating-point issues are solved", {
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 16778071, 16778771
  )
  res <- bed_slop(b, h19, left = 200, right = 200)
  expect_equal(res$start, 16777871)
  expect_equal(res$end, 16778971)
})

## order is different compared to bedtools
test_that("test that old floating-point issues are solved", {
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 160, 170,
    "chr1", 100, 200
  )
  res <- bed_slop(b, h19, both = 0.1, fraction = TRUE)
  expect_equal(res$start, c(90, 159))
  expect_equal(res$end, c(210, 171))
})

test_that("test negative slop on l with strand", {
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 300, 320
  )
  res <- bed_slop(b, tiny.genome, left = -60, right = 60)
  expect_equal(res$start, 360)
  expect_equal(res$end, 380)
})

test_that("test negative slop on l with strand", {
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 300, 320, "a1", 5, "-"
  )
  res <- bed_slop(b, tiny.genome, left = -60, right = 60, strand = TRUE)
  expect_equal(res$start, 240)
  expect_equal(res$end, 260)
})

test_that("test negative slop on r with strand", {
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 300, 320, "a1", 5, "-"
  )

  res <- bed_slop(b, tiny.genome, left = 60, right = -60, strand = TRUE)
  expect_equal(res$start, 360)
  expect_equal(res$end, 380)
})

test_that("test crossover during negative slop", {
  tiny.genome <- tibble::tribble(
    ~ chrom, ~ size,
    "chr1", 1000
  )
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 300, 320, "a1", 5, "-"
  )
  res <- bed_slop(b, tiny.genome, left = -60, right = -60, strand = TRUE)
  expect_equal(res$start, 260)
  expect_equal(res$end, 360)
})

test_that("test crossover during negative slop resulting in 0 length intervals are tossed out", {
  tiny.genome <- tibble::tribble(
    ~ chrom, ~ size,
    "chr1", 1000
  )
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 300, 320, "a1", 5, "-"
  )
  expect_warning(res <- bed_slop(b, tiny.genome, left = -10, right = -10, strand = TRUE))
  expect_equal(nrow(res), 0)
})

test_that("going beyond the end of the chrom", {
  tiny.genome <- tibble::tribble(
    ~ chrom, ~ size,
    "chr1", 1000
  )
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 950, 970, "a1", 5, "-"
  )
  res <- bed_slop(b, tiny.genome, left = 60, right = -60, strand = TRUE, trim = TRUE)
  expect_equal(res$start, 999)
  expect_equal(res$end, 1000)
})

test_that("test edge cases", {
  tiny.genome <- tibble::tribble(
    ~ chrom, ~ size,
    "chr1", 1000
  )
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 50, 60, "a1", 5, "-"
  )
  res <- bed_slop(b, tiny.genome, left = -60, right = 60, strand = TRUE, trim = TRUE)
  expect_equal(res$start, 0)
  expect_equal(res$end, 1)
})

test_that("test edge cases", {
  tiny.genome <- tibble::tribble(
    ~ chrom, ~ size,
    "chr1", 1000
  )
  b <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 50, 60, "a1", 5, "-"
  )
  res <- bed_slop(b, tiny.genome, left = -100, right = 100, strand = TRUE, trim = TRUE)
  expect_equal(res$start, 0)
  expect_equal(res$end, 1)
})
