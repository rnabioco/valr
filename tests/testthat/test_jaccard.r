context("bed_jaccard")

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 10, 20,
  "chr1", 30, 40
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 15, 20
)

test_that("jaccard coeff is calculated correctly", {
  res <- bed_jaccard(x, y)
  expect_equal(res$jaccard, 0.25)
})

test_that("jaccard coeff is calc'd for large data sets", {
  genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

  x <- bed_random(genome, n = 1e5, seed = 10000)
  y <- bed_random(genome, n = 1e5, seed = 20000)

  res <- bed_jaccard(x, y)
  expect_equal(round(res$jaccard, 3), 0.016)
})

test_that("jaccard with grouped inputs are calculated", {
  genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

  x <- bed_random(genome, n = 1e5, seed = 10000)
  y <- bed_random(genome, n = 1e5, seed = 20000)

  res <- bed_jaccard(
    group_by(x, chrom),
    group_by(y, chrom)
  )

  expect_equal(nrow(res), 24)
  expect_true("chrom" %in% names(res))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/jaccard/test-jaccard.sh
test_that("Test symmetry", {
  res <- bed_jaccard(x, y)
  res2 <- bed_jaccard(y, x)
  expect_equal(res$jaccard, res2$jaccard)
})

test_that("Test jaccard with mixed strand files", {
  a <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 10L, 50L, "a1f", 2L, "+",
    "chr1", 20L, 60L, "b1r", 4L, "-",
    "chr1", 25L, 70L, "c1q", 8L, ".",
    "chr1", 30L, 75L, "d1q", 16L, ".",
    "chr1", 40L, 80L, "e1f", 32L, "+",
    "chr1", 45L, 90L, "f1r", 64L, "-",
    "chr2", 10L, 50L, "a2q", 2L, ".",
    "chr2", 20L, 40L, "b2f", 4L, "+",
    "chr2", 25L, 50L, "c2r", 8L, "-",
    "chr2", 30L, 60L, "d2f", 16L, "+",
    "chr2", 35L, 65L, "e2q", 32L, ".",
    "chr2", 39L, 80L, "f2r", 64L, "-"
  )
  b <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 10L, 50L, "2a1r", 2L, "-",
    "chr1", 40L, 70L, "2b1q", 4L, ".",
    "chr1", 60L, 100L, "2c1f", 8L, "+",
    "chr2", 15L, 40L, "2d2f", 16L, "+",
    "chr2", 30L, 100L, "2e2r", 32L, "-"
  )
  res <- bed_jaccard(a, b)
  expect_equal(res$len_i, 145)
  expect_equal(res$len_u, 325)
  expect_equal(round(res$jaccard, 5), round(0.8055556, 5))
  expect_equal(res$n, 2)
})
