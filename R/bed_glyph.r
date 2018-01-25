#' Create example glyphs for valr functions.
#'
#' Used to illustrate the output of valr functions with small examples.
#'
#' @param expr expression to evaluate
#' @param label column name to use for label values. should be present in the
#'   result of the call.
#'
#' @return [ggplot2::ggplot()]
#'
#' @examples
#' x <- trbl_interval(
#'  ~chrom, ~start, ~end,
#'  'chr1', 25,     50,
#'  'chr1', 100,    125
#' )
#'
#' y <- trbl_interval(
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1', 30,     75,  50
#' )
#'
#' bed_glyph(bed_intersect(x, y))
#'
#' x <- trbl_interval(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75,
#'   'chr1', 50,     90,
#'   'chr1', 91,     120
#' )
#'
#' bed_glyph(bed_merge(x))
#'
#' bed_glyph(bed_cluster(x), label = '.id')
#'
#' @export
bed_glyph <- function(expr, label = NULL) {
  expr <- substitute(expr)

  # assign `expr <- quote(bed_intersect(x, y))` at this point to debug
  args_all <- formals(match.fun(expr[[1]]))

  # get required args i.e. those without defaults
  args_req <- names(args_all[sapply(args_all, is.name)])

  # for bed_intersect replace ... with y
  if (expr[[1]] == "bed_intersect") args_req[args_req == "..."] <- "y"

  args_excl <- c("genome", "...")
  args_req <- args_req[!args_req %in% args_excl]

  nargs <- length(args_req)

  # evaluate the expression in the environment context
  env <- parent.frame()
  res <- eval(expr, envir = env)

  # bail if the result is too big
  max_rows <- 100
  if (nrow(res) > max_rows) {
    stop("max_rows exceeded in bed_glyph.", call. = FALSE)
  }

  # get default columns
  cols_default <- c("chrom")
  if ("start" %in% names(res)) cols_default <- c(cols_default, "start")
  if ("end" %in% names(res)) cols_default <- c(cols_default, "end")

  cols_vars <- rlang::syms(cols_default)
  cols_out <- select(res, !!! cols_vars)

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

  # these are the equivalent of the `x` and `y` formals, except are the names
  # of the args in the quoted call.
  expr_vars <- all.vars(expr)

  # this fetches the `x` and `y` rows from the environment
  for (i in 1:nargs) {
    env_i <- get(expr_vars[i], env)
    rows <- mutate(env_i, .facet = expr_vars[i])
    res <- bind_rows(res, rows)
  }

  # assign `.y` values in the result based on clustering
  ys <- group_by(res, .facet)
  ys <- bed_cluster(ys)
  ys <- group_by(ys, .facet, .id)
  ys <- mutate(ys, .y = row_number(.id))
  ys <- ungroup(ys)

  ys <- arrange(ys, .facet, chrom, start)
  res <- arrange(res, .facet, chrom, start)

  res <- mutate(res, .y = ys$.y)

  # make name_result col appear last in the facets
  fct_names <- c(expr_vars, name_result)
  res <- mutate(res, .facet = factor(.facet, levels = fct_names))

  # plot title
  title <- deparse(substitute(expr))

  glyph_plot(res, title, label) + glyph_theme()
}

#' plot for bed_glyph
#' @noRd
glyph_plot <- function(.data, title = NULL, label = NULL) {

  # Colorbrewer 3-class `Greys`
  fill_colors <- c("#f0f0f0", "#bdbdbd", "#636363")

  glyph <- ggplot(.data) +
    geom_rect(
      aes_string(
        xmin = "start", xmax = "end",
        ymin = ".y", ymax = ".y + 0.5",
        fill = ".facet"
      ),
      color = "black", alpha = 0.75
    ) +
    facet_grid(
      .facet ~ .,
      switch = "y",
      scales = "free_y", space = "free_y"
    ) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = NULL, y = NULL)

  if (!is.null(label)) {
    label <- as.name(label)
    aes_label <- aes_(
      x = quote((end - start) / 2 + start),
      y = quote(.y + 0.25),
      label = substitute(label)
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
