x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100, 200,
  "chr1", 150, 250,
  "chr1", 400, 500
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 175, 200,
  "chr1", 175, 225
)

test_that("simple overlap works", {
  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 4)
})

test_that("invert param works", {
  res <- bed_intersect(x, y, invert = TRUE)
  expect_equal(nrow(res), 1)
})

test_that("multiple as", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr1", 100, 200,
    "chr1", 100, 200,
    "chr1", 100, 200,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 175, 200
  )

  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)
})

test_that("multple bs", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 175, 200,
    "chr1", 175, 200,
    "chr1", 175, 200,
    "chr1", 175, 200,
    "chr1", 175, 200
  )

  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 5)
})

test_that("no overlaps returns empty df", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 300, 400
  )
  res <- bed_intersect(x, y)
  expect_true("data.frame" %in% class(res))
  expect_equal(nrow(res), 0)
})

test_that("duplicate intervals are removed (#23)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 500,
    "chr1", 175, 200
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 150, 400,
    "chr1", 151, 401
  )

  res <- bed_intersect(x, y)
  expect_equal(nrow(res), 4)
})

test_that("suffixes disambiguate x/y columns (#28)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1500, ".", ".", "-"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1200, ".", ".", "-"
  )

  res <- bed_intersect(x, y)
  expect_true("start.y" %in% colnames(res))
})

test_that("incorrect `suffix` args throw errors", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score,
    "chr1", 1000, 1500, ".", "."
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score,
    "chr1", 1000, 1200, ".", "."
  )

  expect_error(bed_intersect(x, y, suffix = "TESTING"))
})

test_that("intersections from x bed_tbl with more chroms than y are captured", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr3", 400, 500
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr3", 425, 475
  )

  res <- bed_intersect(x, y)
  expect_true("chr3" %in% res$chrom)
})

test_that("intersections from y bed_tbl with more chroms are captured", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr3", 400, 500
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr3", 425, 475
  )

  res <- bed_intersect(x, y)
  expect_true("chr3" %in% res$chrom)
})

test_that("input x groups are used for comparing intervals issue #108", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "A",
    "chr1", 200, 400, "A",
    "chr1", 300, 500, "A",
    "chr1", 125, 175, "B",
    "chr1", 150, 200, "B"
  )
  x <- arrange(x, chrom, start)
  x <- group_by(x, group)
  res <- bed_intersect(x, x)
  expect_true(all(res$group.x == res$group.y))
})

test_that("tbls grouped by strand are processed", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1500, ".", ".", "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1200, ".", ".", "-"
  )

  res <- bed_intersect(group_by(x, strand), group_by(y, strand))
  expect_equal(nrow(res), 0)

  # flip strands
  res <- bed_intersect(group_by(x, strand), group_by(flip_strands(y), strand))
  expect_equal(nrow(res), 1)
})

test_that("invert = T, and custom suffixes dont result in failed anti_join()", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr3", 500, 600
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 200,
    "chr3", 425, 475
  )

  res <- bed_intersect(x, y, invert = T, suffix = c("a", "b"))
  expect_equal(nrow(res), 1)
})

test_that("multiple y tbl_intervals can be passed to bed_intersect (#220)", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 500,
    "chr2", 200, 400,
    "chr2", 300, 500,
    "chr2", 800, 900
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 350, 430, 300
  )

  z <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 750, 900, 400
  )

  res <- bed_intersect(x, y, z)
  expect_true(all(c("y", "z") %in% res$.source))

  # check that named args can be passed also
  res <- bed_intersect(x, first_file = y, second_file = z)
  expect_true(all(c("first_file", "second_file") %in% res$.source))

  # check that list input is parsed correctly
  res1 <- bed_intersect(x, first_file = y, second_file = z)
  res2 <- bed_intersect(x, list(first_file = y, second_file = z))
  expect_equal(res1, res2)

  res1 <- bed_intersect(x, y, z)
  res2 <- bed_intersect(x, list(y, z))
  expect_equal(res1, res2)
})

test_that("groups are respected when passing multiple y tbl_intervals ", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1500, ".", ".", "+"
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1200, ".", ".", "-"
  )

  z <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 1000, 1200, ".", ".", "+"
  )
  x <- group_by(x, strand)
  y <- group_by(y, strand)
  z <- group_by(z, strand)

  res <- bed_intersect(x, y, z)
  expect_equal(nrow(res), 1)
})

test_that("same intervals are reported with single and multiple intersection", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 500,
    "chr2", 200, 400,
    "chr2", 300, 500,
    "chr2", 800, 900
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 350, 430, 300
  )

  z <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 750, 900, 400
  )
  a <- bed_intersect(x, y)
  b <- bed_intersect(x, z)
  orig <- bind_rows(a, b) |>
    arrange(chrom, start.x, start.y)
  new <- bed_intersect(x, y, z) |>
    arrange(chrom, start.x, start.y) |>
    select(-.source)
  expect_true(all(orig == new))
})

test_that("unmatched groups are included when invert = TRUE", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 500, "A",
    "chr2", 200, 400, "B", # unmatched
    "chr2", 300, 500, "A",
    "chr2", 800, 900, "A"
  ) |> group_by(chrom, group)

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 150, 400, "A",
    "chr1", 500, 550, "A",
    "chr2", 230, 430, "A",
    "chr2", 350, 430, "A"
  ) |> group_by(chrom, group)

  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr2", 200, 400, "B", # unmatched
    "chr2", 800, 900, "A"
  )

  res <- bed_intersect(x, y, invert = TRUE)
  expect_equal(res, pred, ignore_attr = TRUE)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/intersect/test-intersect.sh
test_that("0 length records", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr7", 33059403, 33059403
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group, ~type,
    "chr7", 32599076, 33069221, "NAq", "intron",
    "chr7", 33059336L, 33060883L, "NT5C3A", "intron"
  )
  res <- bed_intersect(x, y)
  expect_equal(res$start.x - res$end.x, c(0, 0))
  expect_equal(res$start.y - res$end.y, c(-470145, -1547))
})


test_that("list input is robustly handled #380", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 100, 500,
    "chr2", 200, 400,
    "chr2", 300, 500,
    "chr2", 800, 900
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 350, 430, 300
  )

  z <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 150, 400, 100,
    "chr1", 500, 550, 100,
    "chr2", 230, 430, 200,
    "chr2", 750, 900, 400
  )

  lst <- list(y, z)

  expect_equal(nrow(bed_intersect(x, y, z)), 11)
  expect_equal(nrow(bed_intersect(x, list(y, z))), 11)
  expect_equal(nrow(bed_intersect(x, lst[1:2])), 11)

  expect_equal(nrow(bed_intersect(x, lst)), 11)
  expect_equal(nrow(bed_intersect(x, lst[[1]], lst[[2]])), 11)
  expect_equal(nrow(bed_intersect(x, lst[1])), 6)
})
