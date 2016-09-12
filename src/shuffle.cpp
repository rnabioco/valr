#include "valr.h"

//[[Rcpp::export]]
DataFrame shuffle_impl(DataFrame df, DataFrame genome, DataFrame interval_bounds,
                       bool within = false, int max_tries = 1000, int seed = 0) {

  // seed for reproducible intervals
  if (seed == 0)
    seed = round(R::runif(0, RAND_MAX)) ;
  
  // seed the generator
  auto generator = ENG(seed) ;
  
  // generated weighted chromosomes to sample from
  CharacterVector chroms_genome = genome["chrom"] ;
  NumericVector sizes_genome = genome["size"] ;
  int nchrom = chroms_genome.size() ;
  
  // calculate weights for chrom distribution
  float mass = sum(sizes_genome) ;
  NumericVector weights = sizes_genome / mass ;
  
  Range chromidx(0, nchrom) ;
  PDIST chrom_dist(chromidx.begin(), chromidx.end(), weights.begin()) ;

  // now deal with query data frame
  CharacterVector chroms_df = df["chrom"] ;
  IntegerVector starts_df = df["start"] ;
  IntegerVector ends_df = df["end"] ;
  
  IntegerVector sizes_df = ends_df - starts_df ;
  
  // number of rows input == rows in output
  int nr = df.size() ;  
  
  CharacterVector chroms_out(nr) ;
  IntegerVector starts_out(nr) ; 
  IntegerVector ends_out(nr) ; 
 
  std::unordered_map<int, intervalTree> interval_trees ;
 
  for (int i = 0; i<nr; ++i) {
    
    // select a chromosome 
    if (within) {
      chroms_out[i] = chroms_df[i] ;
    } else {
      auto chrom_idx = chrom_dist(generator) ;
      chroms_out[i] = chroms_genome[chrom_idx] ;
    }
    
    int rand_start = 0 ; 
    int niter = 0 ;
    bool inbounds = false ;
    intervalVector overlaps ;
    
    while (!inbounds) {
     
      niter++ ;  
   
      tree->findOverlapping(rand_start, rand_end, overlaps) ;
      
      if (!overlaps.empty()) continue ;
        
      if (niter == max_tries) {
        stop('maximum iterations exceeded in bed_shuffle') ;
      }
      
    }  
    
    starts_out[i] = rand_start ;
    ends_out[i] = rand_start + sizes_df[i] ;
  } 
 
  return DataFrame::create( Named("chrom") = chroms_out,
                            Named("start") = starts_out,
                            Named("end") = ends_out) ;
}

/***R
library(dplyr)
library(valr)
genome <- tibble::tribble(
   ~chrom,  ~size,
   "chr1", 500,
   "chr2", 600,
   "chr3", 800
) 

x <- tibble::tribble(
   ~chrom, ~start, ~end,
   "chr1", 100,    300,
   "chr1", 200,    400,
   "chr2", 1,      100,
   "chr2", 200,    400,
   "chr3", 500,    600
) 

x <- bed_shuffle(x, genome)
*/
