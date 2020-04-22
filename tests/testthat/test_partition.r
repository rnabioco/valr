#https://github.com/bedops/bedops/blob/master/applications/bed/bedops/test/TestPlan.xml#L1541
context("bed_partition")

test_that("basic partition works (bedops partition1 test)", {

  x <- tibble::tribble(
    ~chrom,   ~start,   ~end,
    "chr1",   10L,  100L,
    "chr1",   50L,  125L,
    "chr1", 2000L, 2500L,
    "chr0",  250L,  400L,
    "chr1",  250L,  400L,
    "chr1", 2100L, 2125L,
    "chr21",  500L, 1000L,
    "chr0",  100L,  300L,
    "chr1",   50L,  125L,
    "chr1", 2000L, 2500L
  )

  pred <- tibble::tribble(
    ~chrom,   ~start,   ~end,
    "chr0",  100L,  250L,
    "chr0",  250L,  300L,
    "chr0",  300L,  400L,
    "chr1",   10L,   50L,
    "chr1",   50L,  100L,
    "chr1",  100L,  125L,
    "chr1",  250L,  400L,
    "chr1", 2000L, 2100L,
    "chr1", 2100L, 2125L,
    "chr1", 2125L, 2500L,
    "chr21",  500L, 1000L
  )

  res <- bed_partition(x)
  expect_equivalent(res, pred)
})


test_that("extended partition works (bedops partition2 test)", {

x <- tibble::tribble(
  ~chrom,   ~start,   ~end,
  "chr1",   10L,  100L,
  "chr1",   50L,  125L,
  "chr1", 2000L, 2500L,
  "chr3",    1L,    2L,
  "chr0",  100L,  300L,
  "chr1",   50L,  125L,
  "chr1", 2000L, 2500L,
  "chr2",    5L,    7L,
  "chr1",   10L,  100L,
  "chr1",   50L,  125L,
  "chr1",   50L,  125L,
  "chr1", 2000L, 2500L,
  "chr2",    1L,   10L,
  "chr2",    1L,   10L,
  "chr2",    1L,   10L,
  "chr2",    1L,   10L,
  "chr2",    2L,   10L,
  "chr1",    1L,   10L,
  "chr1",    3L,    6L,
  "chr1",    9L,   10L,
  "chr2",    1L,   10L,
  "chr1",    5L,   20L,
  "chr1",   10L,   20L,
  "chr1",   11L,   21L,
  "chr1",   12L,   22L,
  "chr1",   13L,   23L,
  "chr1",   14L,   24L,
  "chr1",   15L,   25L,
  "chr1",   16L,   26L,
  "chr1",   17L,   27L,
  "chr1",   18L,   28L,
  "chr1",   19L,   29L,
  "chr1",   20L,   30L
  )

pred <- tibble::tribble(
  ~chrom,   ~start,   ~end,
  "chr0",  100L,  300L,
  "chr1",    1L,    3L,
  "chr1",    3L,    5L,
  "chr1",    5L,    6L,
  "chr1",    6L,    9L,
  "chr1",    9L,   10L,
  "chr1",   10L,   11L,
  "chr1",   11L,   12L,
  "chr1",   12L,   13L,
  "chr1",   13L,   14L,
  "chr1",   14L,   15L,
  "chr1",   15L,   16L,
  "chr1",   16L,   17L,
  "chr1",   17L,   18L,
  "chr1",   18L,   19L,
  "chr1",   19L,   20L,
  "chr1",   20L,   21L,
  "chr1",   21L,   22L,
  "chr1",   22L,   23L,
  "chr1",   23L,   24L,
  "chr1",   24L,   25L,
  "chr1",   25L,   26L,
  "chr1",   26L,   27L,
  "chr1",   27L,   28L,
  "chr1",   28L,   29L,
  "chr1",   29L,   30L,
  "chr1",   30L,   50L,
  "chr1",   50L,  100L,
  "chr1",  100L,  125L,
  "chr1", 2000L, 2500L,
  "chr2",    1L,    2L,
  "chr2",    2L,    5L,
  "chr2",    5L,    7L,
  "chr2",    7L,   10L,
  "chr3",    1L,    2L
  )

  res <- bed_partition(x)
  expect_equivalent(res, pred)
})


test_that("partition drops non-grouped cols (bedops partition3 test)", {

  x <- tibble::tribble(
   ~chrom,   ~start,   ~end,  ~name, ~score, ~strand, ~seq,
  "chr1", 33657L, 33687L,        "+MA0068.1-Pax4", 8.67655e-06, "+", "TAATGCTATCCCTCCCCCAGCCCCCCACCC",
  "chr1", 33666L, 33686L,       "+MA0073.1-RREB1", 1.97929e-06, "+",           "CCCTCCCCCAGCCCCCCACC",
  "chr1", 33670L, 33690L,       "+MA0073.1-RREB1",  1.0924e-06, "+",           "CCCCCAGCCCCCCACCCACT",
  "chr1", 33672L, 33682L,         "+MA0079.2-SP1", 5.66765e-06, "+",                     "CCCAGCCCCC",
  "chr1", 34375L, 34390L, "+MA0065.2-PPARG::RXRA", 7.21085e-07, "+",                "GGTGGGCAAAGGGCA",
  "chr1", 34377L, 34390L,       "+MA0114.1-HNF4A", 5.44281e-06, "+",                  "TGGGCAAAGGGCA"
  )

  pred <- tibble::tribble(
    ~chrom,   ~start,   ~end,
     "chr1", 33657L, 33666L,
     "chr1", 33666L, 33670L,
     "chr1", 33670L, 33672L,
     "chr1", 33672L, 33682L,
     "chr1", 33682L, 33686L,
     "chr1", 33686L, 33687L,
     "chr1", 33687L, 33690L,
     "chr1", 34375L, 34377L,
     "chr1", 34377L, 34390L
     )

  res <- bed_partition(x)
  expect_equivalent(res, pred)
})


test_that("partition drops non-grouped cols (bedops partition4 test)", {

  x <- tibble::tribble(
    ~chrom,   ~start,   ~end,
  "chr1", 279L, 280L,
  "chr1", 280L, 281L,
  "chr1", 281L, 282L,
  "chr1", 310L, 311L,
  "chr1", 310L, 320L,
  "chr1", 311L, 312L,
  "chr1", 312L, 313L,
  "chr1", 312L, 318L,
  "chr1", 313L, 314L
  )

  pred <- tibble::tribble(
    ~chrom,   ~start,   ~end,
     "chr1", 279L, 280L,
     "chr1", 280L, 281L,
     "chr1", 281L, 282L,
     "chr1", 310L, 311L,
     "chr1", 311L, 312L,
     "chr1", 312L, 313L,
     "chr1", 313L, 314L,
     "chr1", 314L, 318L,
     "chr1", 318L, 320L
     )

  res <- bed_partition(x)
  expect_equivalent(res, pred)
})


test_that("grouping is respected", {

  x <- tibble::tribble(
    ~chrom,   ~start,   ~end,  ~strand,
    "chr1", 33657L, 33687L,"+",
    "chr1", 33666L, 33686L,"+",
    "chr1", 33670L, 33690L,"-",
    "chr1", 33672L, 33682L,"-",
    "chr1", 34375L, 34390L,"+",
    "chr1", 34377L, 34390L,"-"
  )

  x <- group_by(x, strand)
  res <- bed_partition(x)
  expect_true("strand" %in% colnames(res))
  expect_true(all(c("+" %in% res$strand, "-" %in% res$strand)))
  expect_equal(nrow(res), 8)
})

test_that("book-ended intervals are not merged", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100L, 200L,
    "chr1", 200L, 250L
  )

  res <- bed_partition(x)
  expect_equal(res, x)
})

x <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ value, ~id,
  "chr1", 100L, 200L, 1L, "A",
  "chr1", 250L, 500L, 2L, "A",
  "chr2", 250L, 500L, 3L, "A",
  "chr1", 100L, 150L, 10L, "B",
  "chr1", 150L, 250L, 20L, "B",
  "chr1", 140L, 250L, 30L, "B",
  "chr1", 150L, 200L, 40L, "B"
)

test_that("summary functions are executed", {
  res <- bed_partition(x, count = sum(value))
  expect_equal(sum(res$count), 198)
  expect_equal(nrow(res), 6)
})

test_that("summary functions are executed per group", {
  res <- bed_partition(group_by(x, id),
                       count = sum(value, na.rm = T))

  expect_equal(sum(res$count), 196)
  expect_equal(nrow(res), 7)
})

test_that("Tests for multiple columns and operations", {
  res <- bed_partition(x,
                       count = sum(value),
                       max = max(value))
  expect_true(all(c("count", "max") %in% colnames(res)))
  expect_equal(sum(res$count), 198)
  expect_equal(sum(res$max), 115)
})
