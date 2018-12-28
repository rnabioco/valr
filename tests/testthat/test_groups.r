context("dplyr grouping")

genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- read_bed(valr_example("6fields.bed.gz"), n_fields = 6)
y <- x

x_grpd <- group_by(x, strand)
y_grpd <- group_by(y, strand)


test_that("mismatched groups are dropped by two table verbs", {
  res1 <- bed_closest(x, y_grpd)
  res2 <- bed_closest(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 26)

  res1 <- bed_intersect(x, y_grpd)
  res2 <- bed_intersect(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 16)

  res1 <- bed_map(x, y_grpd, .n = n())
  res2 <- bed_map(x_grpd, y, .n = n())
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 10)

  res1 <- bed_subtract(x, y_grpd)
  res2 <- bed_subtract(x_grpd, y)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 0)

  res1 <- bed_window(x, y_grpd, genome)
  res2 <- bed_window(x_grpd, y, genome)
  expect_equal(res1, res2)
  expect_equal(nrow(res1), 16)
})

## dplyr v0.7.9 and earlier style grouped data_frame
df_old <- structure(
    list(
      chrom = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3"),
      start = c(100, 200, 500, 125, 150, 100),
      end = c(200,400, 600, 175, 200, 300),
      group = c("B", "A", "C", "C", "A","A")
    ),
    row.names = c(NA,-6L),
    class = c("grouped_df", "tbl_df",
              "tbl", "data.frame"),
    vars = "group",
    drop = TRUE,
    indices = list(c(1L, 4L, 5L), 0L, 2:3),
    group_sizes = c(3L, 1L, 2L),
    biggest_group_size = 3L,
    labels = structure(
      list(group = c("A", "B", "C")),
      row.names = c(NA,-3L),
      class = "data.frame",
      vars = "group",
      drop = TRUE
    )
  )

## dplyr v0.8.0 style grouped data_frame
df_new <- structure(
    list(
      chrom = c("chr1", "chr1", "chr1", "chr2", "chr2",  "chr3"),
      start = c(100, 200, 500, 125, 150, 100),
      end = c(200,  400, 600, 175, 200, 300),
      group = c("B", "A", "C", "C", "A","A")
    ),
    row.names = c(NA,-6L),
    class = c("grouped_df", "tbl_df",
              "tbl", "data.frame"),
    groups = structure(
      list(group = c("A", "B", "C"),
           .rows = list(c(2L, 5L, 6L), 1L, 3:4)),
      row.names = c(NA,-3L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

test_that("old dataframe groupings (dplyr v. < 0.7.9.900) are tolerated", {

  if (packageVersion("dplyr") >= "0.7.9.9000"){
    expect_warning(bed_intersect(df_old, df_old))
    res <- suppressWarnings(bed_intersect(df_old, df_old))
  } else {
    expect_silent(bed_intersect(df_old, df_old))
    res <- bed_intersect(df_old, df_old)

    #check that input data attributes are not modified
    df_attr_names <-  names(attributes(df_old))
    expect_true(all(c("labels", "indices", "group_sizes") %in% df_attr_names))
    expect_false("groups" %in% df_attr_names)
  }

  res <- suppressWarnings(bed_intersect(df_old, df_old))
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 6)
})


test_that("new style dataframe groupings  (dplyr v. >= 0.7.9.900) are tolerated", {

  # dplyr < 0.7.9.9000 will not recognize the new format
  # groups therefore not respected until regrouped
  if (packageVersion("dplyr") < "0.7.9.9000"){
    res <- bed_intersect(df_new, df_new)
    expect_equal(nrow(res), 10)
    expect_is(res, "data.frame")

    df_new <- group_by(df_new, group)
    res <- bed_intersect(df_new, df_new)
    expect_equal(nrow(res), 6)
    expect_is(res, "data.frame")
  } else {
    res <- bed_intersect(df_new, df_new)
    expect_is(res, "data.frame")
    expect_equal(nrow(res), 6)

    df_new <- ungroup(df_new)
    res <- bed_intersect(df_new, df_new)
    expect_is(res, "data.frame")
    expect_equal(nrow(res), 10)
  }
})


