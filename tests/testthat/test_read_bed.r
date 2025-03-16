bed3_path <- system.file("extdata", "3fields.bed.gz", package = "valr")
bed6_path <- system.file("extdata", "6fields.bed.gz", package = "valr")
bed12_path <- system.file("extdata", "mm9.refGene.bed.gz", package = "valr")
bedgraph_path <- system.file("extdata", "test.bg.gz", package = "valr")
narrowpeak_path <- system.file(
  "extdata",
  "sample.narrowPeak.gz",
  package = "valr"
)
broadpeak_path <- system.file(
  "extdata",
  "sample.broadPeak.gz",
  package = "valr"
)
bigwig_path <- system.file("extdata", "hg19.dnase1.bw", package = "valr")
gtf_path <- system.file("extdata", "hg19.gencode.gtf.gz", package = "valr")

bed3_tbl <- read_bed(bed3_path)

test_that("read BED3", {
  bed3_tbl <- read_bed(bed3_path)
  expect_equal(ncol(bed3_tbl), 3)
  expect_equal(nrow(bed3_tbl), 10)
})


test_that("read BED6", {
  bed6_tbl <- read_bed(bed6_path)
  expect_equal(ncol(bed6_tbl), 6)
  expect_equal(nrow(bed6_tbl), 10)
})

test_that("read BED12", {
  bed12_tbl <- read_bed12(bed12_path)
  expect_equal(ncol(bed12_tbl), 12)
})


test_that("read bedGraph", {
  bedgraph_tbl <- read_bedgraph(bedgraph_path)
  expect_equal(ncol(bedgraph_tbl), 4)
  expect_equal(nrow(bedgraph_tbl), 4)
})

test_that("read narrowPeak", {
  x <- read_narrowpeak(narrowpeak_path)
  expect_equal(ncol(x), 10)
})

test_that("read broadPeak", {
  x <- read_broadpeak(broadpeak_path)
  expect_equal(ncol(x), 9)
})

test_that("read bigwig", {
  x <- read_bigwig(bigwig_path)
  expect_equal(ncol(x), 4)
})


test_that("read gtf", {
  expect_error(read_gtf(gtf_path))
  # expect_equal(ncol(x), 26)
})
