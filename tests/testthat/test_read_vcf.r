context("read_vcf")

v <- system.file('extdata', 'test.vcf.gz', package = 'valr')
x <- read_vcf(v)

test_that('colnames are set in vcf df', {
   test_that('chrom' %in% colnames(x), TRUE)
   test_that('start' %in% colnames(x), TRUE)
   test_that('end' %in% colnames(x), TRUE)
})

test_that('chrom names are set correctly', {
  test_that(all(stringr::str_detect(x$chrom, '^chr')), TRUE)
})

