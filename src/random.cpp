// random.cpp
// 
// generate random intervals on a genome
// 
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Rcpp ;

typedef boost::mt19937 rng_type ;
typedef boost::uniform_int<> dist_type ;
typedef boost::variate_generator<rng_type&, dist_type> gen_type ;

// [[Rcpp::export]]
DataFrame random_impl(DataFrame genome, int length, int n, unsigned int seed = 0) {
 
  CharacterVector chroms = genome["chrom"] ;
  NumericVector sizes = genome["size"] ;
 
  int nchrom = chroms.size() ;
  
  rng_type rng(seed) ;
  dist_type chrom_dist(0, nchrom - 1) ;
  gen_type chrom_rng(rng, chrom_dist) ;
  
  // make and store a RNG for each chrom size
  std::vector< gen_type > size_rngs ;
  
  for (int i=0; i<nchrom; ++i) {
    
    // sub length to avoid off-chrom coordinates
    int size = sizes[i] ;
    dist_type size_dist(1, size - length) ;
    gen_type size_rng(rng, size_dist) ;
    
    size_rngs.push_back(size_rng) ;
  }
  
  CharacterVector rand_chroms(n) ;
  NumericVector rand_starts(n) ;
  
  for (int i=0; i<n; ++i) {
    
     int chrom_idx = chrom_rng() ;
     rand_chroms[i] = chroms[chrom_idx] ;
     
     gen_type size_rng = size_rngs[chrom_idx] ;
     
     int rand_start = size_rng() ;
     rand_starts[i] = rand_start ;
  }
  
  NumericVector rand_ends = rand_starts + length ; 
  
  return DataFrame::create( Named("chrom") = rand_chroms,
                            Named("start") = rand_starts,
                            Named("end") = rand_ends) ;
  
}
