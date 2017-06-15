context('bed_projection')

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1000,
  "chr2", 1000
 )

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 50,    75,
  "chr1", 250,    400,
  "chr1", 500,    600
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,   200,
  "chr1", 200,    800
)

bed_projection(x, y, genome)
bed_projection(x, y, genome, by_chrom = TRUE)

test_that("basic projection test works", {
  #compare output to pbinom( )function
  # 2 of 3 hits, 7 of 10 chance
  # 1- pbinom because pval > .5

  exp <- 1 - pbinom(2, 3, 7/10)
  res <- bed_projection(x, y, genome)
  expect_equal(res$p.value, exp)
})

test_that("projection per chromosome works (by_chrom = TRUE)", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 2000,
    "chr2", 2000
  )

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100,    200,
    "chr1", 250,    400,
    "chr1", 500,    600,
    "chr1", 1000,   2000,
    "chr2", 100,    200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 150,    175,
    "chr1", 525,    575,
    "chr1", 1100,   1200,
    "chr1", 1400,   1600,
    "chr2", 100,    1500
  )
  res <- bed_projection(x, y, genome, by_chrom = TRUE)
  expect_true(all(c("chr1", "chr2") %in% res$chrom))
})

test_that("report significant when intervals are underrepresented, .lower_tail = TRUE", {
  sig <- pbinom(0, 4, .7)

  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 1000,
    "chr2", 1000
  )

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10,    20,
    "chr1", 25,    40,
    "chr1", 50,    60,
    "chr1", 100,   150
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 200,    275,
    "chr1", 525,    575,
    "chr1", 200,    900
  )
  res <- bed_projection(x, y, genome)
  expect_true(res$lower_tail == TRUE)
  expect_true(res$p.value < .5)
  expect_true(res$p.value == sig)
})

test_that("report significant when intervals are overrepresented, .lower_tail = FALSE", {
  sig <- 1 - pbinom(4, 4, .7)
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 1000,
    "chr2", 1000
  )

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 200,    225,
    "chr1", 250,    400,
    "chr1", 500,    600,
    "chr1", 300,   350
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 200,    275,
    "chr1", 525,    575,
    "chr1", 200,    900
  )
  res <- bed_projection(x, y, genome)
  expect_true(res$lower_tail == FALSE)
  expect_true(res$p.value < .5)
  expect_true(res$p.value == sig)
})
