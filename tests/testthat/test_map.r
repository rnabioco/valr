context("bed_map")

test_that("x/y groupings are respected", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ id,
    "chr1", 100, 250, 1,
    "chr2", 250, 500, 2,
    "chr2", 250, 500, 3
  ) %>%
    group_by(id)

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ value, ~ id,
    "chr1", 100, 250, 10, 1,
    "chr1", 150, 250, 20, 2,
    "chr2", 250, 500, 500, 3
  ) %>%
    group_by(id)

  pred <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ id, ~ vals,
    "chr1", 100, 250, 1, 10,
    "chr2", 250, 500, 3, 500,
    "chr2", 250, 500, 2, NA
  )
  res <- bed_map(x, y, vals = sum(value))
  expect_true(all(res == pred, na.rm = TRUE))
})

test_that("values_unique works correctly", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 250
  )

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ value,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20,
    "chr1", 100, 250, 10,
    "chr1", 150, 250, 20
  )

  res <- bed_map(x, y, vals = values_unique(value))
  expect_equal(res$vals, c("10,20"))
})

x <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ id,
  "chr1", 100, 200, 1,
  "chr1", 250, 500, 2,
  "chr2", 250, 500, 3
)

y <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ value,
  "chr1", 100, 150, 10,
  "chr1", 150, 250, 20,
  "chr1", 140, 250, 30,
  "chr1", 150, 200, 40
)

test_that("concat works correctly", {
  res <- bed_map(x, y, vals = concat(value))
  expected <- c("10,30,20,40", NA, NA)
  expect_equal(res$vals, expected)
})

test_that("values works correctly", {
  res <- bed_map(x, y, vals = values(value))
  expected <- c("10,30,20,40", NA, NA)
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
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ value,
    "chr1", 100, 150, 10,
    "chr1", 200, 250, 20
  )

  expected <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ value,
    "chr1", 100, 200, 10
  )
  res <- bed_map(x, y, value = sum(value))
  expect_equivalent(res, expected)
})

test_that("ensure that mapping is calculated with respect to input tbls issue#108", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ group,
    "chr1", 100, 200, "B",
    "chr1", 200, 400, "A",
    "chr1", 500, 600, "C",
    "chr2", 125, 175, "C",
    "chr2", 150, 200, "A",
    "chr3", 100, 300, "A"
  )
  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ group, ~ value,
    "chr1", 100, 199, "A", 10,
    "chr1", 200, 400, "B", 20,
    "chr1", 500, 600, "A", 30,
    "chr2", 125, 175, "C", 40,
    "chr2", 350, 500, "A", 50,
    "chr3", 500, 600, "A", 100
  )

  pred <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ group, ~ total,
    "chr1", 100, 200, "B", NA,
    "chr1", 200, 400, "A", NA,
    "chr1", 500, 600, "C", NA,
    "chr2", 125, 175, "C", 40,
    "chr2", 150, 200, "A", NA,
    "chr3", 100, 300, "A", NA
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group)

  res <- bed_map(x, y, total = sum(value))
  expect_true(all(pred == res, na.rm = T))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/map/test-map.sh
x <- tibble::tribble(
  ~ chrom, ~ start, ~ end,
  "chr1", 0L, 100L,
  "chr1", 100L, 200L,
  "chr2", 0L, 100L,
  "chr2", 100L, 200L,
  "chr3", 0L, 100L,
  "chr3", 100L, 200L
)
y <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ group, ~ value, ~ strand,
  "chr1", 0L, 10L, "a1", 10L, "+",
  "chr1", 10L, 20L, "a2", 5L, "+",
  "chr1", 20L, 30L, "a3", 15L, "+",
  "chr1", 120L, 130L, "a4", 1L, "+",
  "chr3", 0L, 10L, "a5", 1L, "+",
  "chr3", 10L, 20L, "a6", 2L, "+",
  "chr3", 20L, 30L, "a7", 3L, "+",
  "chr3", 120L, 130L, "a8", 4L, "+"
)

## output NA instead of 0, see examples for code to change NA to 0
test_that("test count", {
  res <- bed_map(x, y, vals = n())
  expect_equal(res$vals, c(3, 1, NA, NA, 3, 1))
  res2 <- bed_map(x, y, vals = n()) %>% mutate(vals = ifelse(is.na(vals), 0, vals))
  expect_equal(res2$vals, c(3, 1, 0, 0, 3, 1))
})

# R has no built-in mode function
test_that("test mode", {
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  res <- bed_map(x, y, vals = getmode(value))
  expect_equal(res$vals, c(10, 1, NA, NA, 1, 4))
})

test_that("Test GFF column extraction", {
  z <- tibble::tribble(
    ~ chrom, ~ seqid, ~ type, ~ start, ~ end, ~ score, ~ strand, ~ phase, ~ attributes,
    "chr1", "hg19_ccdsGene", "start_codon", 1L, 9L, 0, "+", ".", "gene_id..CCDS30744.1...transcript_id..CCDS30744.1..",
    "chr1", "hg19_ccdsGene", "CDS", 2L, 11L, 0, "+", "0", "gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";",
    "chr1", "hg19_ccdsGene", "exon", 8L, 20L, 0, "+", ".", "gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";",
    "chr1", "hg19_ccdsGene", "CDS", 9L, 17L, 0, "+", "2", "gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";",
    "chr1", "hg19_ccdsGene", "exon", 40L, 200L, 0, "+", ".", "gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";"
  )

  res <- bed_map(x, z, vals = list(chrom))
  expect_equal(length(res$vals[[1]]), 5)
})

test_that("Tests for multiple columns and operations", {
  res <- bed_map(x, y, count = n(), sum = sum(value))
  expect_equal(res$sum, c(30, 1, NA, NA, 6, 4))
})
