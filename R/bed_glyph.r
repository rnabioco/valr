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
#' bed_glyph2(bed_cluster(x))
#' 
#' @export
bed_glyph <- function(expr, res_name = 'result', label = NULL) {

  excl_params <- c('genome')
  
  expr <- substitute(expr) 
  params <- all.vars(expr) 
  params <- params[! params %in% excl_params]
  nparams <- length(params)
  
  env <- environment()
 
  res <- eval(expr, envir = env)

  if (nparams>1) {
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

  # assign `y` values based on clustering
  res <- bed_cluster(res)
  res <- group_by(res, .facet, .id)
  res <- mutate(res, .y = row_number(.id))
  res <- ungroup(res)
  
  # make res_name col last
  fct_names <- c(params, res_name)
  res <- mutate(res, .facet = factor(.facet, levels = fct_names))
  
  title <- deparse(substitute(expr))
  glyph_plot(res, title, label) + theme_glyph()
}

#' ggplot2 object for glyph
#' @noRd
glyph_plot <- function(res, title = '', label = NULL, colors = c("#fc8d59", "#ffffbf", "#91bfdb")) {
  
  plt <- ggplot(res) + 
    geom_rect(aes_string(xmin = 'start', xmax = 'end',
                         ymin = '.y', ymax = '.y + 0.9',
                         fill = '.facet'),
              color = "black", alpha = 0.75) + 
    facet_grid(.facet ~ ., switch = "y",
               scales = "free_y", space = "free_y") +
    ggtitle(title) +
    scale_fill_manual(values = colors)
  
  if (!is.null(label)) {
    plt <- plt + geom_label(aes_string(x = '(end - start) / 2 + start',
                                       y = '.y + 0.5',
                                       label = 'label')) +
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
