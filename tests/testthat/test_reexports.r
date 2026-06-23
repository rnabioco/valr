test_that("data_frame() is deprecated but still builds a tibble", {
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_warning(
    res <- data_frame(chrom = "chr1", start = 1, end = 10),
    class = "lifecycle_warning_deprecated"
  )
  expect_s3_class(res, "tbl_df")
  expect_equal(res$end, 10)
})

test_that("as_data_frame() is deprecated but still coerces to a tibble", {
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_warning(
    res <- as_data_frame(list(chrom = "chr1", start = 1, end = 10)),
    class = "lifecycle_warning_deprecated"
  )
  expect_s3_class(res, "tbl_df")
})

test_that("read_bigwig() re-export reads example data", {
  skip_if_not_installed("cpp11bigwig")
  res <- read_bigwig(valr_example("hg19.dnase1.bw"))
  expect_true(all(c("chrom", "start", "end") %in% names(res)))
  expect_gt(nrow(res), 0)
})

test_that("read_bigbed() re-export reads example data", {
  skip_if_not_installed("cpp11bigwig")
  res <- read_bigbed(valr_example("test.bb"))
  expect_true(all(c("chrom", "start", "end") %in% names(res)))
  expect_gt(nrow(res), 0)
})
