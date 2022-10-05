x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 75, 125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 50, 100,
  "chr1", 100, 150
)

test_that("reldist calculation is correct", {
  res <- bed_reldist(x, y)
  expect_true(res$.reldist == 0.5)
})

test_that("self reldist is 0", {
  res <- bed_reldist(y, y)
  expect_true(res$.reldist == 0)
})

test_that("detail argument works", {
  res <- bed_reldist(x, y, detail = TRUE)
  expect_true(all(names(res) %in% c("chrom", "start", "end", ".reldist")))
})

test_that("reldist respects groups (#108)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "B",
    "chr1", 200, 400, "A",
    "chr1", 500, 600, "C",
    "chr2", 125, 175, "C",
    "chr2", 150, 200, "A",
    "chr3", 100, 300, "A"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 199, "A",
    "chr1", 200, 400, "B",
    "chr1", 500, 600, "A",
    "chr2", 100, 175, "C",
    "chr2", 350, 500, "A",
    "chr3", 500, 600, "A"
  )

  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 10000,
    "chr2", 10000
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group, chrom)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group, chrom)

  res <- bed_reldist(x, y, detail = TRUE)
  expect_true(nrow(res) == 1)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/reldist/test-reldist.sh
# files not included with package
# a <- read_bed(valr_example("refseq.chr1.exons.bed.gz"), n_fields = 6)
# b <- read_bed(valr_example("aluY.chr1.bed.gz"), n_fields = 6)
# c <- read_bedgraph(valr_example("gerp.chr1.bed.gz"))
#
# test_that("Test intervals that are randomly distributed. The relative distances should equally represented .", {
#   res <- bed_reldist(a, b)
#   d <- tibble::tribble(
#          ~reldist, ~count, ~total, ~fraction,
#                 0,   164L, 43408L,     0.004,
#              0.01,   551L, 43408L,     0.013,
#              0.02,   598L, 43408L,     0.014,
#              0.03,   637L, 43408L,     0.015,
#              0.04,   793L, 43408L,     0.018,
#              0.05,   688L, 43408L,     0.016,
#              0.06,   874L, 43408L,      0.02,
#              0.07,   765L, 43408L,     0.018,
#              0.08,   685L, 43408L,     0.016,
#              0.09,   929L, 43408L,     0.021,
#               0.1,   876L, 43408L,      0.02,
#              0.11,   959L, 43408L,     0.022,
#              0.12,   860L, 43408L,      0.02,
#              0.13,   851L, 43408L,      0.02,
#              0.14,   903L, 43408L,     0.021,
#              0.15,   893L, 43408L,     0.021,
#              0.16,   883L, 43408L,      0.02,
#              0.17,   828L, 43408L,     0.019,
#              0.18,   917L, 43408L,     0.021,
#              0.19,   875L, 43408L,      0.02,
#               0.2,   897L, 43408L,     0.021,
#              0.21,   986L, 43408L,     0.023,
#              0.22,   903L, 43408L,     0.021,
#              0.23,   944L, 43408L,     0.022,
#              0.24,   904L, 43408L,     0.021,
#              0.25,   867L, 43408L,      0.02,
#              0.26,   943L, 43408L,     0.022,
#              0.27,   933L, 43408L,     0.021,
#              0.28,  1132L, 43408L,     0.026,
#              0.29,   881L, 43408L,      0.02,
#               0.3,   851L, 43408L,      0.02,
#              0.31,   963L, 43408L,     0.022,
#              0.32,   950L, 43408L,     0.022,
#              0.33,   965L, 43408L,     0.022,
#              0.34,   907L, 43408L,     0.021,
#              0.35,   884L, 43408L,      0.02,
#              0.36,   965L, 43408L,     0.022,
#              0.37,   944L, 43408L,     0.022,
#              0.38,   911L, 43408L,     0.021,
#              0.39,   939L, 43408L,     0.022,
#               0.4,   921L, 43408L,     0.021,
#              0.41,   950L, 43408L,     0.022,
#              0.42,   935L, 43408L,     0.022,
#              0.43,   919L, 43408L,     0.021,
#              0.44,   915L, 43408L,     0.021,
#              0.45,   934L, 43408L,     0.022,
#              0.46,   843L, 43408L,     0.019,
#              0.47,   850L, 43408L,      0.02,
#              0.48,  1006L, 43408L,     0.023,
#              0.49,   937L, 43408L,     0.022
#          )
#   expect_equal(round(res$.freq, 3), d$fraction)
#   expect_equal(round(res$.counts, 3), d$count)
# })
#
# test_that("Test intervals that are consistently closer to one another than expected.  The distances should be biased towards 0.", {
#   res <- bed_reldist(a, c)
#   d <- tibble::tribble(
#          ~reldist, ~count, ~total, ~fraction,
#                 0, 20629L, 43422L,     0.475,
#              0.01,  2629L, 43422L,     0.061,
#              0.02,  1427L, 43422L,     0.033,
#              0.03,   985L, 43422L,     0.023,
#              0.04,   897L, 43422L,     0.021,
#              0.05,   756L, 43422L,     0.017,
#              0.06,   667L, 43422L,     0.015,
#              0.07,   557L, 43422L,     0.013,
#              0.08,   603L, 43422L,     0.014,
#              0.09,   487L, 43422L,     0.011,
#               0.1,   461L, 43422L,     0.011,
#              0.11,   423L, 43422L,      0.01,
#              0.12,   427L, 43422L,      0.01,
#              0.13,   435L, 43422L,      0.01,
#              0.14,   375L, 43422L,     0.009,
#              0.15,   367L, 43422L,     0.008,
#              0.16,   379L, 43422L,     0.009,
#              0.17,   371L, 43422L,     0.009,
#              0.18,   346L, 43422L,     0.008,
#              0.19,   389L, 43422L,     0.009,
#               0.2,   377L, 43422L,     0.009,
#              0.21,   411L, 43422L,     0.009,
#              0.22,   377L, 43422L,     0.009,
#              0.23,   352L, 43422L,     0.008,
#              0.24,   334L, 43422L,     0.008,
#              0.25,   315L, 43422L,     0.007,
#              0.26,   370L, 43422L,     0.009,
#              0.27,   330L, 43422L,     0.008,
#              0.28,   330L, 43422L,     0.008,
#              0.29,   280L, 43422L,     0.006,
#               0.3,   309L, 43422L,     0.007,
#              0.31,   326L, 43422L,     0.008,
#              0.32,   287L, 43422L,     0.007,
#              0.33,   294L, 43422L,     0.007,
#              0.34,   306L, 43422L,     0.007,
#              0.35,   307L, 43422L,     0.007,
#              0.36,   309L, 43422L,     0.007,
#              0.37,   271L, 43422L,     0.006,
#              0.38,   293L, 43422L,     0.007,
#              0.39,   311L, 43422L,     0.007,
#               0.4,   331L, 43422L,     0.008,
#              0.41,   320L, 43422L,     0.007,
#              0.42,   299L, 43422L,     0.007,
#              0.43,   327L, 43422L,     0.008,
#              0.44,   321L, 43422L,     0.007,
#              0.45,   326L, 43422L,     0.008,
#              0.46,   306L, 43422L,     0.007,
#              0.47,   354L, 43422L,     0.008,
#              0.48,   365L, 43422L,     0.008,
#              0.49,   336L, 43422L,     0.008,
#               0.5,    38L, 43422L,     0.001
#          )
#   expect_equal(round(res$.freq, 3), d$fraction)
#   expect_equal(res$.counts, d$count)
#   expect_equal(res$.total, d$total)
#   expect_equal(sum(res$.counts), sum(d$count))
#
# })
