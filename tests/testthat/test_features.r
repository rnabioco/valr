context("feature functions")

x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

introns <- create_introns(x)
utrs5 <- create_utrs5(x)
utrs3 <- create_utrs3(x)
tss <- create_tss(x)

# helper funs
valid_lengths <- function(x) expect_true(all(x$start < x$end))

test_that("feature lengths are valid", {
  valid_lengths(introns)
  valid_lengths(utrs5)
  valid_lengths(utrs3)
  valid_lengths(tss)
})

test_that("introns and exons don't overlap", {
  # substracting exons from introns should be empty
  expect_equal(nrow(introns %>%
    bed_subtract(x)), 0)
})

test_that("there is 1 feature from each gene", {
  expect_true(all(utrs5 %>%
    group_by(name) %>%
    tally() %>%
    pull(n) == 1))
  expect_true(all(utrs3 %>%
    group_by(name) %>%
    tally() %>%
    pull(n) == 1))
  expect_true(all(tss %>%
    group_by(name) %>%
    tally() %>%
    pull(n) == 1))
})

test_that("TSS are single base features", {
  expect_true(all(tss$end - tss$start == 1))
})
