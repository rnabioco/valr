bed_tbl <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 100, 150, "+",
  "chr1", 200, 250, "+",
  "chr2", 300, 350, "+",
  "chr2", 400, 450, "-",
  "chr3", 500, 550, "-",
  "chr3", 600, 650, "-"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1000,
  "chr2", 2000,
  "chr3", 3000
)

test_that("pos increment works", {
  size <- 100
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(
    all(out$start - bed_tbl$start == size),
    all(out$end - bed_tbl$end == size)
  )
})

test_that("neg increment works", {
  size <- -50
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(
    all(out$start - bed_tbl$start == size),
    all(out$end - bed_tbl$end == size)
  )
})

test_that("starts forced to 0", {
  size <- -120
  out <- bed_shift(bed_tbl, genome, size)
  expect_true(all(out$start >= 0))
})

test_that("end forced to chrom length", {
  size <- 1675
  out <- bed_shift(bed_tbl, genome, size) %>%
    left_join(genome, by = "chrom")
  expect_true(all(out$end <= out$size))
})

test_that("fraction increment works", {
  fraction <- 0.5
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(
    out$start - bed_tbl$start == fraction * interval,
    all(out$end - bed_tbl$end == fraction * interval)
  ))
})

test_that("negative fraction increment works", {
  fraction <- -0.5
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(
    out$start - bed_tbl$start == fraction * interval,
    all(out$end - bed_tbl$end == fraction * interval)
  ))
})

test_that("rounding fraction increment works", {
  fraction <- 0.51234
  interval <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(bed_tbl, genome, fraction = fraction)
  expect_true(all(
    out$start - bed_tbl$start == round(fraction * interval),
    all(out$end - bed_tbl$end == round(fraction * interval))
  ))
})

test_that("shift by strand works", {
  size <- 100
  x <- group_by(bed_tbl, strand)
  out <- bed_shift(x, genome, size)
  expect_true(all(
    ifelse(out$strand == "+",
      out$start - bed_tbl$start == size,
      out$start - bed_tbl$start == -size
    ),
    ifelse(out$strand == "+",
      out$end - bed_tbl$end == size,
      out$end - bed_tbl$end == -size
    )
  ))
})

test_that("shift by strand and fraction works", {
  fraction <- 0.5
  x <- group_by(bed_tbl, strand)
  sizes <- bed_tbl$end - bed_tbl$start
  out <- bed_shift(x, genome, fraction = fraction)
  expect_true(all(
    ifelse(out$strand == "+",
      out$start - bed_tbl$start == sizes * fraction,
      out$start - bed_tbl$start == -sizes * fraction
    ),
    ifelse(out$strand == "+",
      out$end - bed_tbl$end == sizes * fraction,
      out$end - bed_tbl$end == -sizes * fraction
    )
  ))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/shift/test-shift.sh
a <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 100, 200, "a1", 1, "+",
  "chr1", 100, 200, "a2", 2, "-"
)

tiny.genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1000
)

test_that("test going beyond the start of the chrom", {
  out <- bed_shift(a, tiny.genome, size = -300, trim = TRUE)
  expect_true(all(
    out$start == c(0, 0),
    out$end == c(1, 1)
  ))
})

test_that("test going beyond the start of the chrom", {
  out <- bed_shift(a, tiny.genome, size = -200, trim = TRUE)
  expect_true(all(
    out$start == c(0, 0),
    out$end == c(1, 1)
  ))
})

test_that("test going beyond the end of the chrom", {
  out <- bed_shift(a, tiny.genome, size = 1000, trim = TRUE)
  expect_true(all(
    out$start == c(999, 999),
    out$end == c(1000, 1000)
  ))
})

test_that("test shift being larger than a signed int", {
  out <- bed_shift(a, tiny.genome, size = 3000000000, trim = TRUE)
  expect_true(all(
    out$start == c(999, 999),
    out$end == c(1000, 1000)
  ))
})

test_that("test chrom boundaries", {
  tiny2.genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 10
  )

  b <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 5, 10, "cds1", 0, "+"
  )
  out <- bed_shift(b, tiny2.genome, size = 2, trim = TRUE)
  expect_true(all(
    out$start == 7,
    out$end == 10
  ))
})

test_that("test shift huge genome", {
  tiny2.genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 249250621
  )

  b <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 66999638L, 67216822L, "NM_032291", 0L, "+",
    "chr1", 92145899L, 92351836L, "NR_036634", 0L, "-"
  )
  out <- bed_shift(b, tiny2.genome, size = 1000, trim = TRUE)
  expect_true(all(
    out$start == c(67000638, 92146899),
    out$end == c(67217822, 92352836)
  ))
})
