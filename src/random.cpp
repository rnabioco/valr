#include <Rcpp.h>
using namespace Rcpp ;

//[[Rcpp::plugins(cpp11)]]
#include <random>

typedef std::mt19937                       ENG ;
typedef std::uniform_int_distribution<int> DIST;

// [[Rcpp::export]]
DataFrame random_impl(DataFrame genome, int length, int n, int seed = 0) {
 
  CharacterVector chroms = genome["chrom"] ;
  IntegerVector sizes = genome["size"] ;
 
  int nchrom = chroms.size() ;
  
  if (seed == 0)
    seed = round(R::runif(0, RAND_MAX)) ;

  // seed the generator
  auto generator = ENG(seed) ;
 
  DIST chrom_dist(0, nchrom - 1) ;
  
  // make and store a DIST for each chrom size
  std::vector< DIST > size_rngs ;
  
  for (int i=0; i<nchrom; ++i) {
    
    auto size = sizes[i] ;
    // sub length to avoid off-chrom coordinates
    DIST size_dist(1, size - length) ;
    size_rngs.push_back(size_dist) ;
  }
  
  CharacterVector rand_chroms(n) ;
  IntegerVector rand_starts(n) ;
  
  for (int i=0; i<n; ++i) {
     
     auto chrom_idx = chrom_dist(generator) ;
     rand_chroms[i] = chroms[chrom_idx] ;
     
     DIST size_dist = size_rngs[chrom_idx] ;
     
     auto rand_start = size_dist(generator) ;
     rand_starts[i] = rand_start ;
  }
  
  IntegerVector rand_ends = rand_starts + length ; 
  
  return DataFrame::create( Named("chrom") = rand_chroms,
                            Named("start") = rand_starts,
                            Named("end") = rand_ends) ;
  
}

/***R
library(dplyr)
genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 19182237,
  "chr2", 17127713,
  "chr3", 11923987
)

random_impl(genome, length = 1000, n = 10, seed = 0) %>% tibble::as_data_frame()
*/
