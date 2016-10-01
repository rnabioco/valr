#' provides working directory for valr example files
#' @param path path to package files
#' @export
valr_example <- function(path) {
  # https://twitter.com/JennyBryan/status/780150538654527488
  system.file("extdata", path, package = 'valr', mustWork = TRUE)
}


#' convienence function for generating example bed glyphs
#' @param x tbl of intervals passed to .fun
#' @param y optional tbl of intervals to pass to .fun  (default = FALSE)
#' @param .fun valr function to execute on (x and/or y intervals)
#' @param ... other arguments to pass to .fun (x ivls and y ivls passed by default if supplied) 
#' @param x_name name of plotted x intervals (default is "X")
#' @param y_name name of plotted y intervals (default is "Y")
#' @param start_col name of column containing start positions to use for plotting results (default = "start")
#' @param end_col name of columm containing end positions to use for plotting results (default = "end")
#' @param title_ops options to pass to title display in ggplot object (default = "(x, y"))
#' @return \code{ggplot2_obj}
#' @export
bed_glyph <- function(x, y = NULL, 
                      .fun, ...,
                      x_name = "X",
                      y_name = "Y",
                      start_col = "start",
                      end_col = "end",
                      title_ops = "(x, y)") {

  if (!is.null(y)){
    x <- bed_cluster(x)
    x <- mutate(x, .id = as.character(.id))
    x <- group_by(x, .id) 
    x <- mutate(x, bin = row_number(.id),
                title = x_name)
    y <- bed_cluster(y) 
    y <- mutate(y, .id = as.character(.id))
    y <- group_by(y, .id)  
    y <- mutate(y, bin = row_number(.id),
                title = y_name)
    x <- ungroup(x)
    y <- ungroup(y)
    res <- eval(.fun(x, y, ...)) 
    res <- mutate(res, bin = 1, 
                  title = "Result")
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
    x <- mutate(x, .id = as.character(.id))
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
    
    res <- mutate(res, .id = ifelse(".id" %in% colnames(res), 
                                          as.character(.id), 
                                    "1"))
    comb <- bind_rows(x, res) 
    comb <- mutate(comb, title = factor(title,
                                        levels = c(x_name,
                                                   "Result")))
  }
  p_title = paste0(substitute(.fun), title_ops)
  plot_glyph(comb, p_title)
}

#### plotting utility function (not exported)
plot_glyph <- function(p_data, p_title){
  ggplot2::ggplot(p_data) + 
    ggplot2::geom_rect(aes(xmin = start, xmax = end,
                  ymin = bin, ymax = bin + 0.9,
                  fill = title), color = "black", alpha = 0.75) + 
    ggtitle(p_title) +
    facet_grid(title ~ ., switch = "y", scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("#fc8d59", 
                                 "#ffffbf", 
                                 "#91bfdb")) +
    theme_bw() +
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

#' reformat bed tbl to match another tbl
#' 
#' \code{format_bed} returns a tbl whose columns are ordered by another tbl.
#'Columns not found in y tbl are dropped. Missing y columns are added and populated
#' with a dummy entry "."
#' @param x tbl of intervals
#' @param y tbl of intervals 
#' @export
format_bed <- function(x, y) {
  names_x <- names(x)
  names_y <- names(y)
  
  if (any(!names_y %in% names_x)){
    cols_to_add <- setdiff(names_y, names_x)
    n <- ncol(x)
    for (i in seq_along(cols_to_add)){
      x[n + i] <- "."
      colnames(x)[n + i] <- cols_to_add[i]
    }
  }
  x <- select(x, one_of(names_y))

}
