library(dplyr)

#' shift will move each interval in a BED formatted data frame by a user-defined number of bases
#' 
#' @param input_bed data frame in BED format
#' @param genome data frame of chromosome sizes
#' @param increment number of bases to shift intervals, or percent shift if \code{pct=TRUE}
#' @param dna logical indicating shift according to +/- in strand column
#' @param pct logical that defines \code{increment} as a fraction of interval size
#' 
#' @return \code{data.frame}
#' 
#' @seealso \url{http://bedtools.readthedocs.org/en/latest/content/tools/shift.html}
#' 
#' @examples
#' 
#' input_bed <- data.frame(chrom = c("chr1","chr1","chr2","chr2","chr3","chr3"),
#'                          start = c(100,200,300,400,500,600),
#'                          stop = c(150,250,350,450,550,650),
#'                          strand = c('+','+','+','-','-','-'))
#' 
#' genome <- data.frame(chrom = c("chr1","chr2","chr3"),
#'                        size = c(1000,2000,3000))
#' 
#' bed_shift(input_bed,genome,100)
#' bed_shift(input_bed,genome,0.5,pct = TRUE)
#' bed_shift(input_bed,genome,100,dna = TRUE)
#' bed_shift(input_bed,genome,0.5,dna = TRUE,pct = TRUE)
#' 
#' @export

#### define slop function
bed_shift <- function(input_bed, genome, increment, dna=FALSE, pct=FALSE){

  # correct input test
  if (!'data.frame' %in% class(input_bed)){stop('input_bed must be a data frame', call. = FALSE)}
  if (!'data.frame' %in% class(genome)){stop('genome must be a data frame', call. = FALSE)}
  if (!is.numeric(increment)){stop('increment must be numeric', call. = FALSE)}
  if (!is.logical(dna)){stop('dna must be logical', call. = FALSE)}
  if (!is.logical(pct)){stop('pct must be logical', call. = FALSE)}
  if (dna == TRUE && !'strand' %in% colnames(input_bed)){
    stop('To use dna=TRUE, there must be a column named strand in input_bed', call. = FALSE)}
  if (dna == TRUE && any(!'+' %in% input_bed$strand, !'-' %in% input_bed$strand, !length(unique(input_bed$strand)) == 2)){
    stop('strand column must contain only + and -', call. = FALSE)}
  
  # keep original colnames
  original_bed_colnames <- colnames(input_bed)
  original_genome_colnames <- colnames(genome)
  
  # change colnames
  colnames(input_bed)[1:3] <- c("chrom","start","stop")
  colnames(genome) <- c("chrom","chr_len")
  
  # shift invervals
  if (dna == FALSE && pct == FALSE){
    increment <- round(increment)
    input_bed$start <- input_bed$start + increment
    input_bed$stop <- input_bed$stop + increment
  }
  
  # shift by percent of interval size
  if (dna == FALSE && pct == TRUE){
    input_bed <- mutate(input_bed, interval_change = round((input_bed$stop - input_bed$start) * increment) )
    input_bed$start <- input_bed$start + input_bed$interval_change
    input_bed$stop <- input_bed$stop + input_bed$interval_change
    input_bed <- select(input_bed, -interval_change)
  }
  
  # shift by strand
  if (dna == TRUE && pct == FALSE){
    increment <- round(increment)
    input_bed$start[input_bed$strand == '+'] <- input_bed$start[input_bed$strand == '+'] + increment
    input_bed$stop[input_bed$strand == '+'] <- input_bed$stop[input_bed$strand == '+'] + increment
    input_bed$start[input_bed$strand == '-'] <- input_bed$start[input_bed$strand == '-'] - increment
    input_bed$stop[input_bed$strand == '-'] <- input_bed$stop[input_bed$strand == '-'] - increment
  }
  
  # shift by strand and percent
  if (dna == TRUE && pct == TRUE){
    input_bed <- mutate(input_bed, interval_change = round((input_bed$stop - input_bed$start) * increment) )
    
    input_bed$start[input_bed$strand == '+'] <- input_bed$start[input_bed$strand == '+'] + input_bed$interval_change[input_bed$strand == '+']
    input_bed$stop[input_bed$strand == '+'] <- input_bed$stop[input_bed$strand == '+'] + input_bed$interval_change[input_bed$strand == '+']
    input_bed$start[input_bed$strand == '-'] <- input_bed$start[input_bed$strand == '-'] - input_bed$interval_change[input_bed$strand == '-']
    input_bed$stop[input_bed$strand == '-'] <- input_bed$stop[input_bed$strand == '-'] - input_bed$interval_change[input_bed$strand == '-']
    
    input_bed <- select(input_bed, -interval_change)
  }
  
  # if stop < 1, remove row
  if ( any(input_bed$stop < 1) ) {
    input_bed <- input_bed[ -(which(input_bed$stop < 1)) ,]
    warning("Some stop values were < 1 and rows have been removed", call. = FALSE)
  }
  
  # if start < 0, force to 0
  if ( any(input_bed$start < 0) ) {
    input_bed$start[ input_bed$start < 0 ] <- 0
    warning("Some start values were < 0 and have been coerced to 0", call. = FALSE)
  }
  
  # if start >= chr length, remove row
  input_bed <- left_join(input_bed,genome, by = "chrom")
  if ( any(input_bed$start >= input_bed$chr_len)) {
    input_bed <- input_bed[ -(which(input_bed$start >= input_bed$chr_len)) ,]
    warning("Some start values were >= chr length and rows have been removed", call. = FALSE)
  }
  
  # if end past chr length, force to chr length
  if ( any(input_bed$stop > input_bed$chr_len)) {
    input_bed$stop[input_bed$stop > input_bed$chr_len] <- input_bed$chr_len[input_bed$stop > input_bed$chr_len]
    warning("Some stop values were > chrom length and have been coerced to chrom length", call. = FALSE)
  }
  
  # remove chr_len column
  input_bed <- select(input_bed, -chr_len)
  
  # replace original colnames
  colnames(input_bed) <- original_bed_colnames
  
  return(input_bed)
}


#test shift
# bed_shift(input_bed,genome,100)
# bed_shift(input_bed,genome,-120) # with warning that starts were coerced to 0
# bed_shift(input_bed,genome,-220) # with warnings that starts were coerced to 0, and rows were removed
# bed_shift(input_bed,genome,1675) # with warnings that stops were coerced to chrom length, and rows were removed
# bed_shift(input_bed,genome,0.9) # increment is rounded when pct=FALSE
# 
# # test shift by percent
# bed_shift(input_bed,genome,0.5,pct = TRUE)
# bed_shift(input_bed,genome,0.51234,pct = TRUE) # start and stop are rounded when pct=TRUE
# bed_shift(input_bed,genome,2,pct = TRUE)
# bed_shift(input_bed,genome,-5,pct = TRUE) # with warning that stop < 1, and rows removed
# 
# # test shift by strand
# bed_shift(input_bed,genome,100,dna = TRUE)
# bed_shift(input_bed,genome,-220,dna = TRUE) # with warnings that starts were coerced to 0, and rows were removed
# 
# # test shift by strand and percent
# bed_shift(input_bed,genome,0.5,dna = TRUE,pct = TRUE)
# bed_shift(input_bed,genome,-5,dna = TRUE,pct = TRUE)
# 
# # test input checking
# bed_shift(4,genome,100) # all following give input errors
# bed_shift(input_bed,4,100)
# bed_shift(input_bed,genome,'100')
# bed_shift(input_bed,genome,100,dna = 'TRUE')
# bed_shift(input_bed,genome,100,pct = 'TRUE')
# bed_shift(select(input_bed,-strand),genome,100,dna = TRUE) # with warning that strand must be in colnames
# bed_shift(mutate(input_bed, strand = c('+','+','+','-','-','P')),genome,100,dna = TRUE) # strand column must be + -
