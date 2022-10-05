genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- read_bed(valr_example("6fields.bed.gz"), n_fields = 6)
y <- x

x_grpd <- group_by(x, strand)
y_grpd <- group_by(y, strand)


test_that("mismatched groups are dropped by two table verbs", {
  res1 <- bed_closest(x, y_grpd)
  res2 <- bed_closest(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 26)

  res1 <- bed_intersect(x, y_grpd)
  res2 <- bed_intersect(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 16)

  res1 <- bed_map(x, y_grpd, .n = n())
  res2 <- bed_map(x_grpd, y, .n = n())
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 10)

  res1 <- bed_subtract(x, y_grpd)
  res2 <- bed_subtract(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 0)

  res1 <- bed_window(x, y_grpd, genome)
  res2 <- bed_window(x_grpd, y, genome)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 16)
})
