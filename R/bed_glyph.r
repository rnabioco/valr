#' Sample glyph for valr functions
#' 
#' @param expr expression to evaluate
#' @param res_name name of result in output
#' @param label col name for label output
#' 
#' @return \code{ggplot}
#'
#' @examples
#' x <- tibble::tribble(
#'  ~chrom, ~start, ~end,
#'  'chr1', 25,     50,
#'  'chr1', 100,    125
#' )
#'
#' y <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75
#' )
#' 
#' z <- tibble::tribble(
#'   ~chrom, ~start, ~end,
#'   'chr1', 30,     75,
#'   'chr1', 50,     90
#' )
#' 
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 200
#' ) 
#'
#' bed_glyph2(bed_intersect(x, y))
#' bed_glyph2(bed_merge(z))
#' 
#' @export
bed_glyph2 <- function(expr, res_name = 'result', label = NULL) {

  expr <- substitute(expr) 
  params <- all.vars(expr) 
  nparams <- length(params)
  
  env <- environment()
 
  res <- eval(expr, envir = env)

  if (nparams>1 && !'genome' %in% params) {
    res <- select(res, chrom, ends_with(params[1]))
    names_one <- names(get(params[1], env))
    names(res) <- names_one
  }
  
  res <- mutate(res, .facet = res_name)
  
  for (i in 1:nparams) {
    env_i <- get(params[i], env)
    rows <- mutate(env_i, .facet = params[i])
    res <- bind_rows(res, rows)
  } 

  # make res_name col last
  fct_names <- c(params, res_name)
  res <- mutate(res, .facet = factor(.facet, levels = fct_names))
  
  title <- deparse(substitute(expr))
  glyph_plot(res, title, label = label) + theme_glyph()
}

#' ggplot2 object for glyph
#' @noRd
glyph_plot <- function(.data, .title, label, colors = c("#fc8d59", "#ffffbf", "#91bfdb")) {
  
  bin <- 1
  
  plt <- ggplot(.data) + 
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = bin, ymax = bin + 0.9, fill = .facet),
                  color = "black", alpha = 0.75) + 
    facet_grid(.facet ~ ., switch = "y",
               scales = "free_y", space = "free_y") +
    ggtitle(.title) +
    scale_fill_manual(values = colors)
  
  if (!is.null(label)) {
    plt <- plt + geom_label(x = aes((end - start) / 2 + start,
                            y = bin + 0.5,
                            label = label)) +
      xlab(NULL) + ylab(NULL)
  }
  
  plt
}

#' theme for glyphs 
#' @noRd
theme_glyph <- function(base_size = 12, base_family = "Helvetica") {
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
