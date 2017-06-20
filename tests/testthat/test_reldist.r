context("bed_reldist")

x <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    75,       125
)

y <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",    50,       100,
  "chr1",    100,       150
)

test_that("reldist calculation is correct", {
  res <- bed_reldist(x, y)
  expect_true(res$.reldist == 0.5)
})

test_that("self reldist is 0", {
  res <- bed_reldist(y, y)
  expect_true(res$.reldist == 0)
})

test_that("detail argument works", {
  res <- bed_reldist(x, y, detail = TRUE)
  expect_true(all(names(res) %in% c("chrom", "start", "end", ".reldist")))
})

test_that("reldist respects groups (#108)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100,    200,  "B",
    "chr1", 200,    400,  "A",
    "chr1", 500,    600,  "C",
    "chr2", 125,    175,  "C",
    "chr2", 150,    200,  "A",
    "chr3", 100,    300,  "A"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100,    199,  "A",
    "chr1", 200,    400,  "B",
    "chr1", 500,    600,  "A",
    "chr2", 100,    175,  "C",
    "chr2", 350,    500,  "A",
    "chr3", 500,    600,  "A"
  )

  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 10000,
    "chr2", 10000
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group, chrom)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group, chrom)

  res <- bed_reldist(x, y, detail = TRUE)
  expect_true(nrow(res) == 1)

})
