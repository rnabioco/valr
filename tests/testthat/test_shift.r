context('bed_shift')

bed_tbl <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 100, 150, "+",
  "chr1", 200, 250, "+",
  "chr2", 300, 350, "+",
  "chr2", 400, 450, "-",
  "chr3", 500, 550, "-",
  "chr3", 600, 650, "-" )

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1000,
  "chr2", 2000,
  "chr3", 3000 )

test_that("pos increment works", {
  size <- 100
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(all(out$start - bed_tbl$start == size),
              all(out$end - bed_tbl$end == size))
})

test_that("neg increment works", {
  size <- -50
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(all(out$start - bed_tbl$start == size),
              all(out$end - bed_tbl$end == size))
})

test_that("starts forced to 0", {
  size <- -120
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(all(out$start >= 0))
})

test_that("end forced to chrom length", {
  size <- 1675
  out <- bed_shift(bed_tbl, genome, size) %>% left_join(genome, by = "chrom")
  expect_true(all(out$end <= out$size))
})

test_that("fraction increment works", {
  fraction <- 0.5
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(out$start - bed_tbl$start == fraction*interval,
              all(out$end - bed_tbl$end == fraction*interval)))
})

test_that("negative fraction increment works", {
  fraction <- -0.5
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(out$start - bed_tbl$start == fraction*interval,
              all(out$end - bed_tbl$end == fraction*interval)))
})

test_that("rounding fraction increment works", {
  fraction <- 0.51234
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(out$start - bed_tbl$start == round(fraction*interval),
              all(out$end - bed_tbl$end == round(fraction*interval))))
})

test_that("shift by strand works", {
  size <- 100
  x <- group_by(bed_tbl, strand)
  out <- bed_shift(x, genome, size)
  expect_true(all(ifelse(out$strand == '+',
                         out$start - bed_tbl$start == size,
                         out$start - bed_tbl$start == -size),
                  ifelse(out$strand == '+',
                         out$end - bed_tbl$end == size,
                         out$end - bed_tbl$end == -size)))
})

test_that("shift by strand and fraction works", {
  fraction <- 0.5
  x <- group_by(bed_tbl, strand)
  sizes <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(x, genome, fraction = fraction)
  expect_true(all(ifelse(out$strand == '+',
                         out$start - bed_tbl$start == sizes * fraction,
                         out$start - bed_tbl$start == -sizes * fraction),
                  ifelse(out$strand == '+',
                         out$end - bed_tbl$end == sizes * fraction,
                         out$end - bed_tbl$end == -sizes * fraction)))
})
