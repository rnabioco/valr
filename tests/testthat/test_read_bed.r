context("reading BED files")

bed3_path <- system.file('extdata', '3fields.bed.gz', package = 'valr')
bed6_path <- system.file('extdata', '6fields.bed.gz', package = 'valr')
bed12_path <- system.file('extdata', 'mm9.refGene.bed.gz', package = 'valr')
bedgraph_path <- system.file('extdata', 'test.bg.gz', package = 'valr')
narrowpeak_path <- system.file('extdata', 'sample.narrowPeak.gz', package = 'valr')
broadpeak_path <- system.file('extdata', 'sample.broadPeak.gz', package = 'valr')

test_that('read_bed BED3', {
  expect_message(bed3_tbl <- read_bed(bed3_path), "n_fields set automatically to 3")
  expect_equal(ncol(bed3_tbl), 3)
  expect_equal(nrow(bed3_tbl), 10)
})

test_that('read_bed BED6', {
  expect_message(bed6_tbl <- read_bed(bed6_path), "n_fields set automatically to 6")
  expect_equal(ncol(bed6_tbl), 6)
  expect_equal(nrow(bed6_tbl), 10)
})

test_that('read_bed BED12', {
  expect_message(bed12_tbl <- read_bed(bed12_path), "n_fields set automatically to 12")
  expect_equal(ncol(bed12_tbl), 12)
})

test_that('read_bed12', {
  read_bed12_tbl <- read_bed12(bed12_path)
  expect_equal(ncol(read_bed12_tbl), 12)
})

test_that('read bedGraph', {
  bedgraph_tbl <- read_bedgraph(bedgraph_path)
  expect_equal(ncol(bedgraph_tbl), 4)
  expect_equal(nrow(bedgraph_tbl), 4)
})

test_that('read narrowPeak', {
  x <- read_narrowpeak(narrowpeak_path)
  expect_equal(ncol(x), 10)
})

test_that('read broadPeak', {
  x <- read_broadpeak(broadpeak_path)
  expect_equal(ncol(x), 9)
})
