test_that("merge on 1 chrom", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr1", 150, 250,
    "chr1", 200, 350
  )

  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 1)
})

test_that("merge with interval at start", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 50,
    "chr1", 100, 200,
    "chr1", 150, 250
  )

  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 2)
})

test_that("merge with two chroms", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 50,
    "chr1", 25, 75,
    "chr2", 100, 200,
    "chr2", 150, 250
  )

  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 2)
})

test_that("book-ended intervals are merged", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 50,
    "chr1", 50, 100
  )

  res <- bed_merge(bed_df)
  expect_equal(nrow(res), 1)
})

test_that("max_dist is enforced", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 50,
    "chr1", 50, 100
  )

  res <- bed_merge(bed_df, max_dist = 1)
  expect_equal(nrow(res), 1)
})

test_that("max_dist is a positive value", {
  # fmt: skip
  bed_df <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 50,
    "chr1", 50, 100
  )

  expect_error(bed_merge(bed_df, max_dist = -1))
})

test_that("input groups are maintained in the output tbl issue #108", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "A",
    "chr1", 200, 400, "A",
    "chr1", 300, 500, "A",
    "chr1", 125, 175, "C",
    "chr1", 150, 200, "B"
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group)
  res <- bed_merge(x)
  expect_true(all(x$group %in% res$group))
})

test_that("intervals can be merged by strand", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~strand,
    "chr1", 100, 200, "+",
    "chr1", 200, 400, "+",
    "chr1", 300, 500, "+",
    "chr1", 125, 175, "-",
    "chr1", 150, 200, "-"
  )

  x <- group_by(x, strand)
  res <- bed_merge(x)
  expect_equal(nrow(res), 2)
})

test_that("summaries can be computed issue #132", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~value, ~strand,
    "chr1", 1, 50, 1, "+",
    "chr1", 100, 200, 2, "+",
    "chr1", 150, 250, 3, "-",
    "chr2", 1, 25, 4, "+",
    "chr2", 200, 400, 5, "-",
    "chr2", 400, 500, 6, "+",
    "chr2", 450, 550, 7, "+"
  )

  res <- bed_merge(x, .value = sum(value))
  expect_true(all(res$.value != "."))
  expect_true(all(res$.value == c(1, 5, 4, 18)))
})

test_that("multiple summaries can be computed issue #132", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~value, ~strand,
    "chr1", 1, 50, 1, "+",
    "chr1", 100, 200, 2, "+",
    "chr1", 150, 250, 3, "-",
    "chr2", 1, 25, 4, "+",
    "chr2", 200, 400, 5, "-",
    "chr2", 400, 500, 6, "+",
    "chr2", 450, 550, 7, "+"
  )

  res <- bed_merge(x, .value = sum(value), .min = min(value))
  expect_true(all(res$.value != "."))
  expect_true(all(res$.value == c(1, 5, 4, 18)))
  expect_true(all(res$.min == c(1, 2, 4, 5)))
})

test_that("contained intervals are merged issue #176", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 10,
    "chr1", 2, 5,
    "chr1", 7, 9
  )

  res <- bed_merge(x)
  expect_true(nrow(res) == 1)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/merge/test-merge.sh
test_that("Test that precision default is high enough for formatting not to give scientific notation", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand, ~val1, ~val2,
    "chr1", 5333587L, 5344172L, "line1", 0L, "-", 5334680L, 5344172L,
    "chr1", 5481008L, 5484749L, "line2", 0L, "-", 5481796L, 5484749L,
    "chr1", 5481008L, 5484749L, "line3", 0L, "-", 5481796L, 5484749L,
    "chr1", 5481008L, 5484749L, "line4", 0L, "-", 5481796L, 5484749L,
    "chr1", 6763278L, 6766882L, "line5", 0L, "-", 7766544L, 6766882L
  )

  res <- bed_merge(x, .value = sum(val2))
  expect_equal(res$.value, c(5344172, 16454247, 6766882))
})

test_that("Test stranded merge with bedPlus files that have strand", {
  skip_if(packageVersion("readr") <= "1.4.0")

  expect_warning(
    x <- read_bed(valr_example("bug254_e.bed"), skip = 1, lazy = FALSE)
  )
  x <- x |> group_by(strand)
  res <- bed_merge(x, 200) |> arrange(end)
  expect_equal(res$end, c(20000, 25000))
})


test_that("check for off by one errors, related to issue #401 @kcamnairb ", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      10,
    "chr1", 5,      20,
    "chr1", 30,     40
  )
  res <- bed_merge(x, max_dist = 10)
  # fmt: skip
  ex <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      40
  )
  expect_equal(res, ex)

  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      3,
    "chr1", 2,      4,
    "chr1", 5,      10,
    "chr1", 12,     14
  )

  # fmt: skip
  ex <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      4,
    "chr1", 5,      10,
    "chr1", 12,     14
  )
  res <- bed_merge(x, max_dist = 0)
  expect_equal(res, ex)

  # fmt: skip
  ex <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1,      10,
    "chr1", 12,     14
  )
  res <- bed_merge(x, max_dist = 1)
  expect_equal(res, ex)

  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "scaffold_66", 27262, 70396,
    "scaffold_66", 66594, 67647,
    "scaffold_66", 82218, 85280,
    "scaffold_66", 85878, 87553,
    "scaffold_66", 87831, 89885,
    "scaffold_66", 90498, 91996
  )

  # fmt: skip
  ex <- tibble::tribble(
    ~chrom, ~start, ~end,
    "scaffold_66", 27262, 91996
  )
  res <- bed_merge(x, max_dist = 20000)
  expect_equal(res, ex)
})
