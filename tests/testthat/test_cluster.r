context("bed_cluster")

# https://github.com/arq5x/bedtools2/blob/master/test/cluster/test-cluster.sh

x <- tibble::tribble(
  ~chrom, ~start,  ~end,    ~name, ~id, ~strand,
  "chr1", 72017,   884436,  "a",   1,   "+",
  "chr1", 72017,   844113,  "b",   2,   "+",
  "chr1", 939517,  1011278, "c",   3,   "+",
  "chr1", 1142976, 1203168, "d",   4,   "+",
  "chr1", 1153667, 1298845, "e",   5,   "-",
  "chr1", 1153667, 1219633, "f",   6,   "+",
  "chr1", 1155173, 1200334, "g",   7,   "-",
  "chr1", 1229798, 1500664, "h",   8,   "-",
  "chr1", 1297735, 1357056, "i",   9,   "+",
  "chr1", 1844181, 1931789, "j",   10,  "-"
)

test_that("basic cluster works", {
  res <- bed_cluster(x)
  # test number of groups in output
  expect_equal(length(unique(res$.id)), 4)
})

test_that("stranded cluster works", {
  res <- bed_cluster(group_by(x, strand))
  # test number of groups in output
  expect_equal(length(unique(res$.id)), 6)
})

x <- tibble::tribble(
  ~chrom, ~start,  ~end,    ~name, ~id, ~strand,
  "chr1", 72017,   884436,  "a",   1,   "+",
  "chr1", 72017,   844113,  "b",   2,   "+",
  "chr1", 939517,  1011278, "c",   3,   "+",
  "chr2", 940000,  990000,  "d",   4,   "-"
)

test_that("cluster ids are not repeated per group issue #171", {
  res <- bed_cluster(x)
  # test that groups have unique ids
  chr1_ids <- filter(res, chrom == "chr1") %>%
    select(.id) %>% unique() %>%  unlist()
  chr2_ids <- filter(res, chrom == "chr2") %>%
    select(.id) %>% unique() %>%  unlist()
  shared_ids <- intersect(chr1_ids, chr2_ids)
  expect_equal(length(shared_ids), 0)
})
