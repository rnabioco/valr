context("astyle")

test_that("source code formatting", {
  skip("run `astyle` test manually")
  skip_on_cran()
  skip_on_os("windows")
  skip_on_travis()
  skip_on_appveyor()

  expect_warning(astyle("--dry-run"), NA)
})
