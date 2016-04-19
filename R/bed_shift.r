#' shift will move each interval in a BED formatted data frame by a user-defined number of bases
#' 
#' @param input_bed data frame in BED format
#' @param genome data frame of chromosome sizes
#' @param increment number of bases to shift intervals to the right, or negative increment to shift to the left, or percent shift if \code{pct=TRUE}
#' @param dna logical indicating shift according to +/- in strand column
#' @param pct logical that defines \code{increment} as a fraction of interval size
#' 
#' @return \code{data.frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/shift.html}
#' 
#' @examples
#' 
#' input_bed <- tibble(
#'    ~chrom, ~start, ~end, ~strand,
#'    "chr1", 100, 150, "+",
#'    "chr1", 200, 250, "+",
#'    "chr2", 300, 350, "+",
#'    "chr2", 400, 450, "-",
#'    "chr3", 500, 550, "-",
#'    "chr3", 600, 650, "-" )
#' 
#' genome <- tibble(
#'    ~chrom, ~size,
#'    "chr1", 1000,
#'    "chr2", 2000,
#'    "chr3", 3000 )
#' 
#' bed_shift(input_bed,genome,100)
#' bed_shift(input_bed,genome,0.5,pct = TRUE)
#' bed_shift(input_bed,genome,100,dna = TRUE)
#' bed_shift(input_bed,genome,0.5,dna = TRUE,pct = TRUE)
#' 
#' @export

# define slop function
bed_shift <- function(input_bed, genome, increment, dna=FALSE, pct=FALSE){

  # correct input test
  if (!'data.frame' %in% class(input_bed)){stop('input_bed must be a data frame', call. = FALSE)}
  if (!'data.frame' %in% class(genome)){stop('genome must be a data frame', call. = FALSE)}
  if (!is.numeric(increment)){stop('increment must be numeric', call. = FALSE)}
  if (!is.logical(dna)){stop('dna must be logical', call. = FALSE)}
  if (!is.logical(pct)){stop('pct must be logical', call. = FALSE)}
  if (dna && !'strand' %in% colnames(input_bed)){
    stop('To use dna=TRUE, there must be a column named strand in input_bed', call. = FALSE)}
  if (dna && any(!'+' %in% input_bed$strand, !'-' %in% input_bed$strand, !length(unique(input_bed$strand)) == 2)){
    stop('strand column must contain only + and -', call. = FALSE)}
  if ( !all(colnames(input_bed)[1:3] == c("chrom","start","end")) ){
    stop('first three input_bed colnames must be "chrom" "start" "end"', call. = FALSE)}
  if ( !all(colnames(genome) == c("chrom","size")) ){
    stop('genome colnames must be "chrom" "size"', call. = FALSE)}
  
  # work on out_bed
  out_bed <- input_bed
  
  # shift invervals
  if (!dna && !pct){
    out_bed <- out_bed %>%
      mutate(start = start + round(increment),
             end = end + round(increment))
  }
  
  # shift by percent of interval size
  if (!dna && pct){
    out_bed <- out_bed %>%
      mutate(interval_change = end - start) %>%
      mutate(start = start + round(interval_change * increment),
             end = end + round(interval_change * increment)) %>%
      select(-interval_change)
  }
  
  # shift by strand
  if (dna && !pct){
    out_bed <- out_bed %>%
      mutate(start = ifelse(strand == '+',
                            start + round(increment),
                            start - round(increment)),
             end = ifelse(strand == '+',
                          end + round(increment),
                          end - round(increment)))
  }
  
  # shift by strand and percent
  if (dna && pct){
    out_bed <- out_bed %>%
      mutate(interval_change = end - start) %>%
      mutate(start = ifelse(strand == "+",
                            start + round(interval_change * increment),
                            start - round(interval_change * increment)),
             end = ifelse(strand == "+",
                          end + round(interval_change * increment),
                          end - round(interval_change * increment))) %>%
      select(-interval_change)
  }
  
  # genome size checking
  out_bed <- out_bed %>%
    left_join(genome, by = "chrom") %>%
    filter(start < size, end >= 1) %>%
    mutate(start = ifelse(start < 0, 0, start),
           end = ifelse(end > size, size, end)) %>%
    select(-size)
  
  return(out_bed)
}
