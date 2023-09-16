# https://github.com/arq5x/bedtools2/blob/master/test/cluster/test-cluster.sh

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~id, ~strand,
  "chr1", 72017, 884436, "a", 1, "+",
  "chr1", 72017, 844113, "b", 2, "+",
  "chr1", 939517, 1011278, "c", 3, "+",
  "chr1", 1142976, 1203168, "d", 4, "+",
  "chr1", 1153667, 1298845, "e", 5, "-",
  "chr1", 1153667, 1219633, "f", 6, "+",
  "chr1", 1155173, 1200334, "g", 7, "-",
  "chr1", 1229798, 1500664, "h", 8, "-",
  "chr1", 1297735, 1357056, "i", 9, "+",
  "chr1", 1844181, 1931789, "j", 10, "-"
)

test_that("basic cluster works", {
  res <- bed_cluster(x)
  # test number of groups in output
  expect_equal(length(unique(res$.id)), 4)
  expect_equal(res$.id, c(1, 1, 2, 3, 3, 3, 3, 3, 3, 4))
})

test_that("stranded cluster works", {
  res <- bed_cluster(group_by(x, strand))
  # test number of groups in output
  expect_equal(length(unique(res$.id)), 6)
  expect_equal(res$.id, c(1, 1, 2, 3, 3, 4, 4, 4, 5, 6))
})

x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~id, ~strand,
  "chr1", 72017, 884436, "a", 1, "+",
  "chr1", 72017, 844113, "b", 2, "+",
  "chr1", 939517, 1011278, "c", 3, "+",
  "chr2", 940000, 990000, "d", 4, "-"
)

test_that("cluster ids are not repeated per group issue #171", {
  res <- bed_cluster(x)
  # test that groups have unique ids
  chr1_ids <- filter(res, chrom == "chr1") |>
    select(.id) |>
    unique() |>
    unlist()
  chr2_ids <- filter(res, chrom == "chr2") |>
    select(.id) |>
    unique() |>
    unlist()
  shared_ids <- intersect(chr1_ids, chr2_ids)
  expect_equal(length(shared_ids), 0)
})


test_that("guard against max_dist argument preventing clustering first interval in contig issue #388", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "a", 1, 10,
    "a", 20, 50,
    "b", 20, 50,
    "c", 100, 100
  )

  res <- bed_cluster(x, max_dist = 0)
  expect_equal(res$.id, 1L:4L)

  res <- bed_cluster(x, max_dist = 100)
  expect_equal(res$.id, c(1, 1, 2, 3))

  res <- bed_cluster(x, max_dist = 10)
  expect_equal(res$.id, c(1, 1, 2, 3))

  res <- bed_cluster(x, max_dist = 9)
  expect_equal(res$.id, 1L:4L)
})

test_that("check for off by one errors, related to issue #401 @kcamnairb ", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      10,
    "chr1", 5,      20,
    "chr1", 30,     40
  )
  res <- bed_cluster(x, max_dist = 10)
  expect_equal(res$.id, c(1L, 1L, 1L))

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      3,
    "chr1", 2,      4,
    "chr1", 5,      10,
    "chr1", 12,     14
  )
  res <- bed_cluster(x, max_dist = 0)
  expect_equal(res$.id, c(1L, 1L, 2L, 3L))

  res <- bed_cluster(x, max_dist = 1)
  expect_equal(res$.id, c(1L, 1L, 1L, 2L))
})

test_that("check for additional errors, related to issue #401 @kcamnairb ", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "scaffold_66", 27262, 70396,
    "scaffold_66", 66594, 67647,
    "scaffold_66", 82218, 85280,
    "scaffold_66", 85878, 87553,
    "scaffold_66", 87831, 89885,
    "scaffold_66", 90498, 91996
  )
  res <- bed_cluster(x, max_dist = 20000)
  expect_true(all(res$.id == 1))

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 10,
    "chr1", 1, 11,
    "chr1", 1, 9,
    "chr1", 1, 9,
    "chr1", 3, 4,
    "chr1", 3, 12,
    "chr1", 10, 14,
    "chr1", 100, 200,
    "chr2", 1, 10,
    "chr2", 1, 11,
    "chr2", 1, 9,
    "chr2", 1, 9,
    "chr2", 3, 4,
    "chr2", 3, 12,
    "chr2", 10, 14,
    "chr2", 100, 200
  )

  res <- bed_cluster(x, max_dist = 0)
  expect_true(max(res$.id) == 4)

  res <- bed_cluster(x, max_dist = 100)
  expect_equal(res$.id, c(rep(1, 8), rep(2, 8)))

  res <- bed_cluster(x, max_dist = -3)
  expect_true(max(res$.id) == 6)
})
