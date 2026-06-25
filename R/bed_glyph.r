#' Create example glyphs for valr functions.
#'
#' Used to illustrate the output of valr functions with small examples.
#'
#' @param expr expression to evaluate
#' @param label column name to use for label values. should be present in the
#'   result of the call.
#' @param max_rows maximum number of rows in the evaluated result that can be
#'   plotted. Calls producing more rows than this raise an error.
#'
#' @return [ggplot2::ggplot()]
#'
#' @examples
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 25,     50,
#'   "chr1", 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~value,
#'   "chr1", 30, 75, 50
#' )
#'
#' bed_glyph(bed_intersect(x, y))
#'
#' x <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   "chr1", 30,     75,
#'   "chr1", 50,     90,
#'   "chr1", 91,     120
#' )
#'
#' bed_glyph(bed_merge(x))
#'
#' bed_glyph(bed_cluster(x), label = ".id")
#'
#' @export
bed_glyph <- function(expr, label = NULL, max_rows = 100L) {
  expr <- substitute(expr)

  # evaluate the expression in the calling environment
  env <- parent.frame()
  res <- eval(expr, envir = env)

  # bail if the result is too big
  if (nrow(res) > max_rows) {
    cli::cli_abort(
      "Result has {nrow(res)} row{?s}, exceeding {.arg max_rows} ({max_rows}).
       Use a smaller example or raise {.arg max_rows}."
    )
  }

  # input interval tbls are the positional (unnamed) arguments of the call that
  # are bare symbols resolving to an interval data frame in `env`, e.g. `x` and
  # `y`. Named arguments (`min_overlap = 0L`) and non-interval positionals (a
  # `genome` tbl, which has no start/end) are ignored. This avoids per-function
  # special-casing and works for namespaced calls like `valr::bed_merge(x)`.
  call_args <- as.list(expr)[-1]
  arg_names <- names(call_args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(call_args))
  }
  positional <- call_args[arg_names == ""]
  is_input <- vapply(
    positional,
    function(a) {
      if (!is.symbol(a) || !exists(as.character(a), envir = env)) {
        return(FALSE)
      }
      obj <- get(as.character(a), envir = env)
      is.data.frame(obj) && all(c("start", "end") %in% names(obj))
    },
    logical(1)
  )
  # unique() collapses repeated args (e.g. `bed_intersect(x, x)`) to one facet
  input_vars <- unique(vapply(positional[is_input], as.character, character(1)))

  # get default columns
  cols_default <- c("chrom")
  if ("start" %in% names(res)) {
    cols_default <- c(cols_default, "start")
  }
  if ("end" %in% names(res)) {
    cols_default <- c(cols_default, "end")
  }

  cols_out <- select(res, all_of(cols_default))

  # get cols that are now suffixed in the result. This is a reasonable default
  # for bed_intersect and functions that call bed_intersect.
  suffix_default <- stringr::fixed(".x")
  cols_out <- bind_cols(cols_out, select(res, ends_with(suffix_default)))

  # get any named columns from the expr
  expr_names <- names(expr)
  expr_names <- expr_names[expr_names != ""]
  expr_names <- intersect(expr_names, names(res))

  if (length(expr_names) > 0) {
    cols_out <- bind_cols(cols_out, select(res, starts_with(expr_names)))
  }

  # get dot cols from result, e.g. `.overlap`
  dot_fixed <- stringr::fixed(".")
  cols_out <- bind_cols(cols_out, select(res, starts_with(dot_fixed)))

  # strip suffixes from names, assumes suffixes are dot-character, e.g. `.x`
  names_strip <- stringr::str_replace(names(cols_out), "\\.[:alnum:]$", "")
  names(cols_out) <- names_strip

  res <- cols_out
  name_result <- "result"
  res <- mutate(res, .facet = name_result)

  # add the input rows, faceted by their argument name
  for (var in input_vars) {
    rows <- mutate(as_tibble(get(var, envir = env)), .facet = var)
    res <- bind_rows(res, rows)
  }

  # assign `.y` (stacking position) from per-facet clustering. Join back by a
  # stable row id rather than relying on two arrange() calls producing identical
  # orders; this also preserves the result's own `.id` column (e.g. from
  # bed_cluster) for use as a `label`.
  res <- mutate(res, .row = row_number())
  ys <- group_by(res, .data[[".facet"]])
  ys <- bed_cluster(ys)
  ys <- group_by(ys, .data[[".facet"]], .data[[".id"]])
  ys <- mutate(ys, .y = row_number())
  ys <- ungroup(ys)
  ys <- select(ys, all_of(c(".row", ".y")))
  res <- left_join(res, ys, by = ".row")
  res <- select(res, -all_of(".row"))

  # make name_result col appear last in the facets
  fct_names <- c(input_vars, name_result)
  res <- mutate(res, .facet = factor(.data[[".facet"]], levels = fct_names))

  # plot title
  title <- deparse(expr)

  glyph_plot(res, title, label) + glyph_theme()
}

#' plot for bed_glyph
#' @noRd
glyph_plot <- function(.data, title = NULL, label = NULL) {
  if (!is.null(label) && !label %in% names(.data)) {
    cli::cli_abort("{.arg label} ({.val {label}}) is not a column in the result.")
  }

  # Colorbrewer 3-class `Greys`
  fill_colors <- c("#f0f0f0", "#bdbdbd", "#636363")

  glyph <- ggplot(.data) +
    geom_rect(
      aes(
        xmin = .data[["start"]],
        xmax = .data[["end"]],
        ymin = .data[[".y"]],
        ymax = .data[[".y"]] + 0.5,
        fill = .data[[".facet"]]
      ),
      color = "black",
      alpha = 0.75
    ) +
    facet_grid(
      rows = vars(.data[[".facet"]]),
      switch = "y",
      scales = "free_y",
      space = "free_y"
    ) +
    # constant (not proportional) y padding so every `.y` unit maps to the same
    # pixel height across facets: all interval glyphs render the same size and
    # the figure height tracks the total row count.
    scale_y_continuous(expand = expansion(add = 0.25)) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = NULL, y = NULL)

  if (!is.null(label)) {
    aes_label <- aes(
      x = (.data[["end"]] - .data[["start"]]) / 2 + .data[["start"]],
      y = .data[[".y"]] + 0.25,
      label = .data[[label]]
    )
    glyph <- glyph + geom_label(aes_label, na.rm = TRUE)
  }

  glyph
}

#' theme for bed_glyph
#' @noRd
glyph_theme <- function(base_size = 12, base_family = "Helvetica") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank()
    )
}
