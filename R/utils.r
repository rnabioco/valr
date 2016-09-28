#' provides working directory for valr example files
#' 
#' @export
valr_example <- function(path) {
  # https://twitter.com/JennyBryan/status/780150538654527488
  system.file("extdata", path, package = 'valr', mustWork = TRUE)
}


#' convienence function for generating example bed glyphs
#'
#' @param expr valr expression 
#' 
#' @return \code{ggplot2_obj}
#'
#' @examples 
#' 
#' bed_glyph(bed_intersect(x, y))
#' 
#' @export
bed_glyph <- function(expr,
                      start_col = "start",
                      end_col = "end") {

  args <- deparse(substitute(expr))
 
  if (!is.null(y)){
    x <- bed_cluster(x)
    x <- group_by(x, .id) 
    x <- mutate(x, bin = row_number(.id),
                title = x_name)
    y <- bed_cluster(y) 
    y <- group_by(y, .id)  
    y <- mutate(y, bin = row_number(.id),
                title = y_name)
    x <- ungroup(x)
    y <- ungroup(y)
    res <- eval(.fun(x, y, ...)) 
    res <- mutate(res, bin = 1, 
                  label = "Result")
    if ("start.x" %in% colnames(res)){
      res <- mutate_(res,
                     start = start_col,
                     end = end_col)
    }
    comb <- bind_rows(x, y, res) 
    comb <- mutate(comb, title = factor(title,
                                        levels = c(x_name,
                                                   y_name,
                                                   "Result")))
  } else if (is.null(y)){
    x <- bed_cluster(x)
    x <- group_by(x, .id) 
    x <- mutate(x, bin = row_number(.id),
                title = x_name)
    x <- ungroup(x)
    res <- eval(.fun(x, ...)) 
    
    
    res <- mutate(res, bin = 1, 
                  title = "Result")
    
    if ("start.x" %in% colnames(res)){
      res <- mutate_(res,
                     start = start_col,
                     end = end_col)
    }
    
    
    comb <- bind_rows(x, res) 
    comb <- mutate(comb, title = factor(title,
                                        levels = c(x_name,
                                                   "Result")))
  }
  title = paste0(substitute(expr))
  glyph_plot(comb, title) + glyph_theme() + glyph_label()
}

#' ggplot2 theme for glyphs
#' 
glyph_theme <- ggplot2::theme(
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
)

#' plotting utility function
#' 
glyph_plot <- function(data, title_plot){
   
  colors_default <- c('#fc8d59', '#ffffbf', '#91bfdb')
  
  ggplot2::ggplot(data) + 
    ggplot2::geom_rect(aes(xmin = start, xmax = end,
                           ymin = bin, ymax = bin + 0.9,
                           fill = title),
                       color = "black", alpha = 0.75) + 
    ggtitle(title_plot) +
    facet_grid(label ~ ., switch = "y",
               scales = "free_y", space = "free_y") +
    scale_fill_manual(values = colors_default) +
    theme_bw() + xlab('') + ylab('')

}

#' labels for glyph plot
#' 
glyph_label <- function(data) {
  ggplot2::geom_label(data, aes((end - start) / 2 + start,
                                bin + 0.5, label = .id)) + 
    xlab("") + ylab("")
}
