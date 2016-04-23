context('bed_shift')

input_bed <- tibble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 100, 150, "+",
  "chr1", 200, 250, "+",
  "chr2", 300, 350, "+",
  "chr2", 400, 450, "-",
  "chr3", 500, 550, "-",
  "chr3", 600, 650, "-" )

genome <- tibble(
  ~chrom, ~size,
  "chr1", 1000,
  "chr2", 2000,
  "chr3", 3000 )

test_that("pos increment works", {
  shift <- 100
  out <- bed_shift(input_bed, genome, shift)
  expect_true(all(out$start - input_bed$start == shift),
              all(out$end - input_bed$end == shift))
})

test_that("neg increment works", {
  shift <- -100
  out <- bed_shift(input_bed, genome, shift)
  expect_true(all(out$start - input_bed$start == shift),
              all(out$end - input_bed$end == shift))
})

test_that("starts forced to 0", {
  shift <- -120
  out <- bed_shift(input_bed, genome, shift)
  expect_true(all(out$start >= 0))
})

test_that("end forced to chrom length", {
  shift <- 1675
  out <- bed_shift(input_bed, genome, shift) %>% left_join(genome, by = "chrom")
  expect_true(all(out$end <= out$size))
})

test_that("increment is rounded", {
  shift <- 0.9
  out <- bed_shift(input_bed, genome, shift)
  expect_true(all(out$start - input_bed$start == round(shift)),
              all(out$end - input_bed$end == round(shift)))
})

test_that("percent increment works", {
  shift <- 0.5
  interval <- input_bed$end - input_bed$start
  out <- bed_shift(input_bed, genome, shift, pct = TRUE)
  expect_true(all(out$start - input_bed$start == shift*interval,
              all(out$end - input_bed$end == shift*interval)))
})

test_that("negative percent increment works", {
  shift <- -0.5
  interval <- input_bed$end - input_bed$start
  out <- bed_shift(input_bed, genome, shift, pct = TRUE)
  expect_true(all(out$start - input_bed$start == shift*interval,
              all(out$end - input_bed$end == shift*interval)))
})

test_that("rounding percent increment works", {
  shift <- 0.51234
  interval <- input_bed$end - input_bed$start
  out <- bed_shift(input_bed, genome, shift, pct = TRUE)
  expect_true(all(out$start - input_bed$start == round(shift*interval),
              all(out$end - input_bed$end == round(shift*interval))))
})

test_that("shift by strand works", {
  shift <- 100
  out <- bed_shift(input_bed, genome, shift, dna = TRUE)
  expect_true(all(ifelse(out$strand == '+',
                         out$start - input_bed$start == shift,
                         out$start - input_bed$start == -shift),
                  ifelse(out$strand == '+',
                         out$end - input_bed$end == shift,
                         out$end - input_bed$end == -shift)))
})

test_that("shift by strand and percent works", {
  shift <- 0.5
  interval <- input_bed$end - input_bed$start
  out <- bed_shift(input_bed, genome, shift, dna = TRUE, pct = TRUE)
  expect_true(all(ifelse(out$strand == '+',
                         out$start - input_bed$start == shift*interval,
                         out$start - input_bed$start == -shift*interval),
                  ifelse(out$strand == '+',
                         out$end - input_bed$end == shift*interval,
                         out$end - input_bed$end == -shift*interval)))
})
