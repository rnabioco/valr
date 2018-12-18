context("factors")

genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
seed <- 1010486
x <- bed_random(genome, n = 1e4, seed = seed)
x$strand <- rep(c("+", "_"), 5e3)

x_facs <- x
x_facs$chrom <- factor(x_facs$chrom)
x_facs$strand <- factor(x_facs$strand)

x_grpd <- group_by(x, strand)
x_facs_grpd <- group_by(x_facs, strand)

test_that("factor types for groups are handled the same as character", {
  res_x <- bed_intersect(x_grpd, x_grpd)
  expect_warning(res_xfacs <- bed_intersect(x_facs_grpd, x_facs_grpd))
  expect_true(all(res_xfacs == res_x))
})


test_that("mixing factor and character vectors for grouping works", {
  res_x <- bed_intersect(x_grpd, x_grpd)
  expect_warning(res_mixed <- bed_intersect(x_grpd, x_facs_grpd))
  expect_true(all(res_x == res_mixed))
})


test_that("factors with no entries are handled ", {
  x_empty_groups <- x_facs_grpd %>%
    filter(strand == "+", chrom == "chr1") %>%
    group_by(strand)

  expect_warning(res_x <- bed_intersect(x_facs_grpd, x_empty_groups))
  expect_true(all(res_x$chrom == "chr1"))
  expect_true(all(res_x$strand.x == "+" && res_x$strand.y == "+"))

})

test_that("complex, raw, and other types are not supported", {
  tmp <- x
  tmp$compex <- 1 + 2i
  expect_error(bed_intersect(tmp, tmp))

  tmp <- x
  tmp$raw <- as.raw(42)
  expect_error(bed_intersect(tmp, tmp))

})

test_that("list columns are supported", {
  x$lst_col <- list(1:10)
  res <- bed_intersect(x, x)
  expect_equal(nrow(res), 10076)
})
