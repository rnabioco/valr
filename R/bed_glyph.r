#' Create example glyphs for valr functions.
#' 
#' Used to illustrate the output of valr functions with small input tbls.
#' 
#' @param expr expression to evaluate
#' @param label colname in output to use for label values
#' @param res_name name of result in output
#'   
#' @return a \code{ggplot} object
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
#' bed_glyph(bed_intersect(x, y))
#' 
#' x <- tibble::tribble(
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
  env <- parent.frame()
  res <- eval(expr, envir = env)

  # need to figure out what columns will be in the result.
  
  # get default columns
  cols_default <- c('chrom')
  if ('start' %in% names(res)) cols_default <- c(cols_default, 'start')
  if ('end' %in% names(res)) cols_default <- c(cols_default, 'end')
  
  out_cols <- select_(res, .dots = cols_default)
 
  # get `x` that are now suffixed in the result. While is is possible for 
  # functions that take a suffix argument (bed_intersect), it is not possible for funcs
  # like bed_map that do not take a suffix and call intersect internally. However because `.x`
  # is the intersect default, it should usually work.
  suffix_default <- '.x'
  out_cols <- bind_cols(out_cols, select(res, ends_with(stringr::fixed(suffix_default))))
 
  # get any named columns from the expr
  expr_names <- names(expr)
  expr_names <- expr_names[expr_names != '']
  expr_names <- intersect(expr_names, names(res))
  
  if (!purrr::is_empty(expr_names)) out_cols <- bind_cols(out_cols, select(res, starts_with(expr_names)))
 
  # get dot cols from result e.g. `.overlap`
  out_cols <- bind_cols(out_cols, select(res, starts_with(stringr::fixed('.'))))
 
  # strip suffixes from names, assumes suffixes are dot-character, `.x`
  names_strip <- stringr::str_replace(names(out_cols), '\\.[:alnum:]$', '') 
  names(out_cols) <- names_strip

  res <- out_cols
  res <- mutate(res, .facet = res_name)

  # these are the equivalent of the `x` and `y` formals, except are the names
  # of the args in the quoted call. 
  expr_vars <- all.vars(expr)
  
  # this fetches the `x` and `y` rows from the environment 
  for (i in 1:nargs) {
    env_i <- get(args_req[i], env)
    rows <- mutate(env_i, .facet = expr_vars[i])
    res <- bind_rows(res, rows)
  } 

  # assign `.y` values based on clustering
  ys <- bed_cluster(res)
  ys <- group_by(ys, .facet, .id)
  ys <- mutate(ys, .y = row_number(.id))
  ys <- ungroup(ys)
  
  res <- mutate(res, .y = ys$.y)
  
  # make res_name col last
  fct_names <- c(expr_vars, res_name)
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
    glyph <- glyph + geom_label(aes_label, na.rm = TRUE)
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
