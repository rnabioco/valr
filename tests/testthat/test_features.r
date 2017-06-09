context("feature functions")

x <- read_bed12(valr_example('mm9.refGene.bed.gz'))

introns <- create_introns(x)
utrs5 <- create_utrs5(x)
utrs3 <- create_utrs3(x)

# helper funs
valid_lengths <- function(x) expect_true(all(x$start < x$end))

test_that("introns are valid", {
  valid_lengths(introns)
})

test_that("5´ UTRs are valid", {
  valid_lengths(utrs5)
})

test_that("3´ UTRs are valid", {
  valid_lengths(utrs3)
})
