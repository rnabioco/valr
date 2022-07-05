# https://github.com/arq5x/bedtools2/blob/master/test/closest/test-closest.sh
context("bed_closest")


test_that("1bp closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1", 20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)
  expect_true(all(c(0, -1, 1) %in% res$.dist))
})

test_that("reciprocal test of 1bp closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10,
    "chr1", 19, 20,
    "chr1", 20, 21
  )
  res <- bed_closest(y, x)
  expect_equal(nrow(res), 3)
  expect_true(all(c(0, -1, 1) %in% res$.dist))
})

test_that("0bp apart closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1", 20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 3)
  expect_true(all(c(0, -1, 1) %in% res$.dist))
})

test_that("reciprocal of 0bp apart closer, check for off-by-one errors", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1", 20, 21
  )
  res <- bed_closest(y, x)
  res2 <- bed_closest(x, y)
  expect_equal(nrow(res), 3)
  expect_equal(nrow(res), 3)
  expect_true(all(c(0, -1, 1) %in% res$.dist))
})

test_that("check that first left interval at index 0 is not lost", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1)
})

test_that("check that first right interval at index 0 is not lost", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 20, 21
  )
  res <- bed_closest(x, y)
  expect_equal(nrow(res), 1)
})

test_that("check that strand closest works (strand = TRUE)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10, "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 90, 120, "b", 1, "-"
  )

  res <- bed_closest(group_by(x, strand), group_by(y, strand))
  expect_equal(nrow(res), 0)
})

test_that("check that same strand is reported (strand = TRUE", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 80, 100, "q1", 1, "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 5, 15, "d1.1", 1, "+",
    "chr1", 20, 60, "d1.2", 2, "-",
    "chr1", 200, 220, "d1.3", 3, "-"
  )

  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.distance,
    "chr1", 80, 100, "q1", 1, "+", 5, 15, "d1.1", 1, "+", 0, -66
  )

  res <- bed_closest(group_by(x, strand), group_by(y, strand))
  expect_true(all(pred == res))
})

test_that("check that different strand is reported (strand_opp = TRUE", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 80, 100, "q1", 1, "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 5, 15, "d1.1", 1, "+",
    "chr1", 20, 60, "d1.2", 2, "-",
    "chr1", 200, 220, "d1.3", 3, "-"
  )

  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "q1", 1, "+", 20, 60, "d1.2", 2, "+", 0, -21
  )

  res <- bed_closest(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_true(all(pred == res))
})

test_that("check that reciprocal strand closest works (strand_opp = TRUE) ", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 100, 200, "a", 10, "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 80, 90, "b", 1, "-"
  )

  res <- bed_closest(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_equal(nrow(res), 1)
})

test_that("overlapping intervals are removed (overlap = F)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 9, 10,
    "chr1", 19, 21,
    "chr1", 20, 21
  )

  res <- bed_closest(x, y, overlap = FALSE)
  expect_true(res[2, "start.y"] != 19)
})

test_that("duplicate intervals are not reported", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr1", 150, 200,
    "chr1", 550, 580,
    "chr2", 7000, 8500
  )
  res <- bed_closest(x, y)
  expect_false(any(duplicated(res)))
})

test_that("all overlapping features are reported", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr1", 150, 200,
    "chr1", 50, 100,
    "chr1", 200, 300
  )
  exp <- tibble::tribble(
    ~chrom, ~start.x, ~start.y,
    "chr1", 100, 200
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 4)
})

test_that("test reporting of first overlapping feature and
           overlap = F excludes overlapping intervals", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 101,
    "chr1", 200, 201,
    "chr1", 300, 301,
    "chr1", 100000, 100010,
    "chr1", 100020, 100040,
    "chr2", 1, 10,
    "chr2", 20, 30
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 101,
    "chr1", 150, 201,
    "chr1", 175, 375
  )
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~start.y, ~end.y, ~.dist,
    "chr1", 100, 101, 150, 201, 50,
    "chr1", 200, 201, 100, 101, -100,
    "chr1", 300, 301, 150, 201, -100,
    "chr1", 100000, 100010, 175, 375, -99626,
    "chr1", 100020, 100040, 175, 375, -99646
  )
  res <- bed_closest(x, y, overlap = F)
  expect_true(all(pred == res))
})

### test distance reporting conditions ###

### tbls to test
d_q1 <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 80, 100, "d_q1.1", 5, "+"
)

d_q2 <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 80, 100, "d_q2.1", 5, "-"
)

d_d1F <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 40, 60, "d1F.1", 10, "+"
)

d_d1R <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 40, 60, "d1R.1", 10, "-"
)

d_d2F <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 140, 160, "d2F.1", 10, "+"
)

d_d2R <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 140, 160, "d2R.1", 10, "-"
)

test_that("default distance reporting works for forward hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q1.1", 5, "+", 40, 60, "d1F.1", 10, "+", 0, -21
  )
  res <- bed_closest(d_q1, d_d1F)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for reverse hit on left, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q1.1", 5, "+", 40, 60, "d1R.1", 10, "-", 0, -21
  )
  res <- bed_closest(d_q1, d_d1R)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for forward hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q2.1", 5, "-", 40, 60, "d1F.1", 10, "+", 0, -21
  )
  res <- bed_closest(d_q2, d_d1F)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for reverse hit on left, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q2.1", 5, "-", 40, 60, "d1R.1", 10, "-", 0, -21
  )
  res <- bed_closest(d_q2, d_d1R)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for forward hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q1.1", 5, "+", 140, 160, "d2F.1", 10, "+", 0, 41
  )
  res <- bed_closest(d_q1, d_d2F)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for reverse hit on right, forward query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q1.1", 5, "+", 140, 160, "d2R.1", 10, "-", 0, 41
  )
  res <- bed_closest(d_q1, d_d2R)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for forward hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q2.1", 5, "-", 140, 160, "d2F.1", 10, "+", 0, 41
  )
  res <- bed_closest(d_q2, d_d2F)
  expect_true(all(pred == res))
})

test_that("default distance reporting works for reverse hit on right, reverse query", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 80, 100, "d_q2.1", 5, "-", 140, 160, "d2R.1", 10, "-", 0, 41
  )
  res <- bed_closest(d_q2, d_d2R)
  expect_true(all(pred == res))
})

### additional tbls for tests ###
a2 <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 10, 20, "a1", 1, "-"
)

b2 <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 8, 9, "b1", 1, "+",
  "chr1", 21, 22, "b2", 1, "-"
)

test_that("Make sure non-overlapping ties are reported ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 10, 20, "a1", 1, "-", 8, 9, "b1", 1, "+", 0, -2,
    "chr1", 10, 20, "a1", 1, "-", 21, 22, "b2", 1, "-", 0, 2
  )
  res <- bed_closest(a2, b2)
  expect_true(all(pred == res))
})

test_that("Make sure non-overlapping ties are reported with strand = T ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 10, 20, "a1", 1, "-", 21, 22, "b2", 1, "-", 0, 2
  )
  res <- bed_closest(group_by(a2, strand), group_by(b2, strand))
  expect_true(all(pred == res))
})

test_that("Make sure non-overlapping ties are reported with strand_opp = T ", {
  pred <- tibble::tribble(
    ~chrom, ~start.x, ~end.x, ~name.x, ~score.x, ~strand.x, ~start.y, ~end.y, ~name.y, ~score.y, ~strand.y, ~.overlap, ~.dist,
    "chr1", 10, 20, "a1", 1, "-", 8, 9, "b1", 1, "-", 0, -2
  )
  res <- bed_closest(group_by(a2, strand), group_by(flip_strands(b2), strand))
  expect_true(all(pred == res))
})

test_that("Make sure that closest intervals are captured when intervals span multiple interval tree nodes issue #105", {
  # when the y tbl has >= 64 intervals two nodes of the interval tree will be generated
  snps <- read_bed(valr_example("hg19.snps147.chr22.bed.gz"), n_fields = 6, n_max = 10)
  genes_one_node <- read_bed(valr_example("genes.hg19.chr22.bed.gz"), n_fields = 6, n_max = 63)
  genes_two_nodes <- read_bed(valr_example("genes.hg19.chr22.bed.gz"), n_fields = 6, n_max = 64)

  res_expt_one_node <- bed_closest(snps, genes_one_node)
  res_expt_two_nodes <- bed_closest(snps, genes_two_nodes)
  # adding one extra interval should not result in doubling the reported intervals
  expect_false(nrow(res_expt_two_nodes) >= 2 * nrow(res_expt_one_node))
})

test_that("test that a max of two duplicated x ivls are returned, assuming non-overlapping, and non-duplicate y ivls #105", {
  # ed
  snps <- read_bed(valr_example("hg19.snps147.chr22.bed.gz"), n_fields = 6, n_max = 10)
  genes <- read_bed(valr_example("genes.hg19.chr22.bed.gz"), n_fields = 6, n_max = 64)
  # make sure there are no repeated y ivls (otherwise more than 2 x ivls should be reported)
  genes <- group_by(genes, chrom, start, end)
  genes <- mutate(genes, ivl_count = n())
  genes <- filter(genes, ivl_count == 1)
  genes <- select(genes, -ivl_count)
  genes <- group_by(genes, chrom)

  res <- bed_closest(snps, genes, overlap = FALSE)
  res <- group_by(res, chrom, start.x, end.x)
  res <- summarize(res, n = n())
  # there should not be more than 2 possible closest ivls.
  expect_true(all(res$n <= 2))
})


test_that("ensure that subtraction is done with respect to input tbls issue#108", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "A",
    "chr1", 200, 400, "A",
    "chr1", 300, 500, "A",
    "chr1", 125, 175, "C",
    "chr1", 150, 200, "B"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "A",
    "chr1", 200, 400, "B",
    "chr1", 300, 500, "A",
    "chr1", 125, 175, "C",
    "chr2", 150, 200, "A"
  )
  x_grouped <- arrange(x, chrom, start) %>%
    group_by(group, chrom)
  y_grouped <- arrange(y, chrom, start) %>%
    group_by(group, chrom)
  res <- bed_closest(x_grouped, y_grouped)
  expect_true(all(res$group.x == res$group.y))
})

# from https://github.com/arq5x/bedtools2/blob/master/test/closest/test-closest.sh
test_that("test closest forcing -s yet no matching strands on chrom", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 100, 200, "a", 10, "+"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 90, 120, "b", 1, "-"
  )
  res <- bed_closest(group_by(x, strand), group_by(y, strand))
  expect_true(nrow(res) == 0)
})

test_that("test closest forcing -S with only an opp strands on chrom", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 100, 200, "a", 10, "+"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 90, 120, "b", 1, "-"
  )
  res <- bed_closest(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_true(nrow(res) == 1)
})

test_that("Make sure non-overlapping ties are reported", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 10, 20, "a1", 1, "-"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 8, 9, "b1", 1, "+",
    "chr1", 21, 22, "b2", 1, "-"
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 2)
})

test_that("Make sure non-overlapping ties are reported, with strand option", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 10, 20, "a1", 1, "-"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 8, 9, "b1", 1, "+",
    "chr1", 21, 22, "b2", 1, "-"
  )
  res <- bed_closest(group_by(x, strand), group_by(y, strand))
  expect_true(nrow(res) == 1)
})

test_that("Make sure non-overlapping ties are reported, with strand-oppo option", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 10, 20, "a1", 1, "-"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 8, 9, "b1", 1, "+",
    "chr1", 21, 22, "b2", 1, "-"
  )
  res <- bed_closest(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_true(nrow(res) == 1)
})

test_that("check ties, single db", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 10, 20, "a1", 1, "-"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~score, ~strand,
    "chr1", 8, 9, "b1", 1, "+",
    "chr1", 21, 22, "b2", 1, "-"
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 2)
})

test_that("check reporting of adjacent intervals issue #348", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 10, 20
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 8, 9,
    "chr1", 9, 10,
    "chr1", 20, 21,
    "chr1", 21, 22
  )
  res <- bed_closest(x, y)
  expect_true(nrow(res) == 2)
})
