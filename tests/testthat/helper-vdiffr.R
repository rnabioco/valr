# from: https://raw.githubusercontent.com/tidyverse/ggplot2/759c63c2fd9e00ba3322c1b74b227f63c98d2e06/tests/testthat/helper-vdiffr.R

# By default, if vdiffr is not installed, all visual tests are skipped unless
# VDIFFR_RUN_TESTS is explicitly set to "true", which should be the case only on
# a GitHub Actions CI runner with stable version of R.

if (
  requireNamespace("vdiffr", quietly = TRUE) &&
    utils::packageVersion("testthat") >= "3.0.3"
) {
  expect_doppelganger <- vdiffr::expect_doppelganger
} else {
  # If vdiffr is not available and visual tests are explicitly required, raise error.
  if (identical(Sys.getenv("VDIFFR_RUN_TESTS"), "true")) {
    abort("vdiffr is not installed")
  }

  # Otherwise, assign a dummy function
  expect_doppelganger <- function(...) skip("vdiffr is not installed.")
}
