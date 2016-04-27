#include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::depends(BH)]]
#include <boost/random.hpp>

// [[Rcpp::export]]
DataFrame random_impl(DataFrame genome, int length, int n, unsigned int seed = 0) {
 
  CharacterVector chroms = genome["chrom"] ;
  NumericVector sizes = genome["size"] ;
 
  int nchrom = chroms.size() ;
  
  if (seed == 0)
    seed = rand() ;
  
  boost::mt19937 rng(seed) ;
  boost::uniform_int<> chrom_dist(0, nchrom - 1) ;
  boost::variate_generator< boost::mt19937&, boost::uniform_int<> > chrom_rng(rng, chrom_dist) ;
  
  // make and store a RNG for each chrom size
  std::vector< boost::variate_generator< boost::mt19937&, boost::uniform_int<> > > size_rngs ;
  
  for (int i=0; i<nchrom; ++i) {
    
    int size = sizes[i] ;
    // sub length to avoid off-chrom coordinates
    boost::uniform_int<> size_dist(1, size - length) ;
    boost::variate_generator< boost::mt19937&, boost::uniform_int<> > size_rng(rng, size_dist) ;
    
    size_rngs.push_back(size_rng) ;
  }
  
  CharacterVector rand_chroms(n) ;
  NumericVector rand_starts(n) ;
  
  for (int i=0; i<n; ++i) {
    
     int chrom_idx = chrom_rng() ;
     rand_chroms[i] = chroms[chrom_idx] ;
     
     boost::variate_generator< boost::mt19937&, boost::uniform_int<> > size_rng = size_rngs[chrom_idx] ;
     
     int rand_start = size_rng() ;
     rand_starts[i] = rand_start ;
  }
  
  NumericVector rand_ends = rand_starts + length ; 
  
  return DataFrame::create( Named("chrom") = rand_chroms,
                            Named("start") = rand_starts,
                            Named("end") = rand_ends) ;
  
}

/*** R
genome <- tibble::frame_data(
  ~chrom,  ~size,
  "chr1",  10000000,
  "chr2",  50000000,
  "chr3",  60000000,
  "chrX",  5000000
)

# random intervals (unsorted)
bed_random(genome)
*/
