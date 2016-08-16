#include "valr.h"

//[[Rcpp::export]]
DataFrame shuffle_impl(DataFrame df, DataFrame incl, DataFrame excl, DataFrame genome,
                       int max_tries) {

  genome_map_t chrom_sizes = makeChromSizes(genome) ;
 
  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;
  
  IntegerVector sizes = ends - starts ; 
  
  CharacterVector chroms = df["chrom"] ;
  
  std::vector<std::string> chroms_out ;
  std::vector<int> starts_out ; 
  std::vector<int> ends_out ; 

  return DataFrame::create( Named("chrom") = chroms_out,
                            Named("start") = starts_out,
                            Named("end") = ends_out) ;
}

/***R
library(dplyr)
library(valr)
genome <- tibble::frame_data(
   ~chrom,  ~size,
   "chr1", 500,
   "chr2", 600,
   "chr3", 800
) 

x <- tibble::frame_data(
   ~chrom, ~start, ~end,
   "chr1", 100,    300,
   "chr1", 200,    400,
   "chr2", 1,      100,
   "chr2", 200,    400,
   "chr3", 500,    600
) 

x <- bed_shuffle(x, genome)
*/
