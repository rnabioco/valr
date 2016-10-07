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
#'   ~chrom, ~start, ~end, ~value,
#'   'chr1', 30,     75,  50
#' )
#' 
#' z <- tibble::tribble(
#'   ~chrom, ~start, ~end, ~label,
#'   'chr1', 30,     75,   '',
#'   'chr1', 50,     90,   '',
#'   'chr1', 91,     120,  'book-ended'
#' )
#' 
#' genome <- tibble::tribble(
#'   ~chrom, ~size,
#'   "chr1", 200
#' ) 
#'
#' bed_glyph(bed_intersect(x, y))
#' bed_glyph(bed_merge(z))
#' bed_glyph(bed_cluster(x))
#' 
#' @export
bed_glyph <- function(expr, label = NULL, res_name = 'result') {

  expr <- substitute(expr)
  
  args_all <- formals(match.fun(expr[[1]]))
  
  # get required args i.e. those without defaults 
  args_req <- names(args_all[sapply(args_all, is.name)])
  
  # remove ellipsis and excl_args
  args_excl <- c('genome')
  args_req <- args_req[args_req != '...']
  args_req <- args_req[!args_req %in% args_excl]
  
  nargs <- length(args_req)
 
  # evaluate the expression in the environment context
  env <- environment()
  res <- eval(expr, envir = env)

  # need to deal with a few different cases:
  # 1. if there is a suffix argument, get the 
  if (nargs>1) {
    if ('suffix' %in% names(args_all)) {
      # get the suffix for the `x` param
      suffix_x <- as.character(args_all$suffix)[2]
      res <- select(res, chrom, ends_with(suffix_x))
    }
    names_x <- names(get(args_req[1], env))
    names(res) <- names_x
  }
  
  res <- mutate(res, .facet = res_name)
  
  for (i in 1:nargs) {
    env_i <- get(args_req[i], env)
    rows <- mutate(env_i, .facet = args_req[i])
    res <- bind_rows(res, rows)
  } 

  # assign `.y` values based on clustering
  res <- bed_cluster(res)
  res <- group_by(res, .facet, .id)
  res <- mutate(res, .y = row_number(.id))
  res <- ungroup(res)
  
  # make res_name col last
  fct_names <- c(args_req, res_name)
  res <- mutate(res, .facet = factor(.facet, levels = fct_names))

  # plotting ------------------------------------------------------- 

  # plot title 
  title <- deparse(substitute(expr))
 
  fill_colors <- c("#fc8d59", "#ffffbf", "#91bfdb")
 
  glyph <- ggplot(res) +
    geom_rect(aes_string(xmin = 'start', xmax = 'end',
                         ymin = '.y', ymax = '.y + 0.9',
                         fill = '.facet'),
              color = "black", alpha = 0.75) + 
    facet_grid(.facet ~ ., switch = "y",
               scales = "free_y", space = "free_y") +
    ggtitle(title) +
    scale_fill_manual(values = fill_colors) +
    xlab(NULL) + ylab(NULL)
  
  if (!is.null(label)) {
    label <- as.name(label)
    aes_label <- aes_(x = quote((end - start) / 2 + start),
                      y = quote(.y + 0.5),
                      label = substitute(label))
    glyph <- glyph + geom_label(aes_label)
  }

  glyph + theme_glyph()
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
