v <- system.file("extdata", "test.vcf.gz", package = "valr")
x <- read_vcf(v)

test_that("colnames are set in vcf df", {
  expect_true("chrom" %in% colnames(x))
  expect_true("start" %in% colnames(x))
  expect_true("end" %in% colnames(x))
})

test_that("chrom names are set correctly", {
  expect_true(all(stringr::str_detect(x$chrom, "^chr")))
})
