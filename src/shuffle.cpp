#include "valr.h"

typedef std::unordered_map<int, intervalTree> chrom_tree_t ;

//[[Rcpp::export]]
DataFrame shuffle_impl(DataFrame df, DataFrame incl, int max_tries = 1000, int seed = 0) {

  // seed for reproducible intervals
  if (seed == 0)
    seed = round(R::runif(0, RAND_MAX)) ;
  
  CharacterVector df_chroms = df["chrom"] ;
  IntegerVector df_starts = df["start"] ;
  IntegerVector df_ends = df["end"] ;
  
  IntegerVector df_sizes = df_ends - df_starts ;
    
  // seed the generator
  auto generator = ENG(seed) ;
 
  int nr = df.nrows() ; 
  CharacterVector chroms_out(nr) ;
  IntegerVector starts_out(nr) ; 
  IntegerVector ends_out(nr) ; 
 
  chrom_tree_t interval_trees ;
 
  for (int i = 0; i<nr; ++i) {
    
    // select a chromosome 
    chroms_out[i] = df_chroms[i] ;
    auto chrom_idx = chrom_dist(generator) ;
    chroms_out[i] = chroms_genome[chrom_idx] ;
    
    int rand_start = 0 ; 
    int niter = 0 ;
    bool inbounds = false ;
    intervalVector overlaps ;
    
    while (!inbounds) {
     
      niter++ ;  
      
      int rand_start = chrom_rng() ;
      int rand_end = df_sizes[i] ;
   
      tree->findContained(rand_start, rand_end, overlaps) ;
     
      if (overlaps.empty()) {
        // didn't find an overlap
        continue ;
      } else if (niter == max_tries) {
        // tried too many times to find an overlap
        stop('maximum iterations exceeded in bed_shuffle') ;
      }
     
      // keep the interval 
      inbounds = true ;
    }  
    
    starts_out[i] = rand_start ;
    ends_out[i] = rand_end ;
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
