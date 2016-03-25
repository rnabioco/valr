// bed_random
// 
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Rcpp ;

//' create random intervals in a genome
//' 
//' @param genome tbl of chrom sizes
//' @param length length of intervals
//' @param n number of intervals to generate
//' @param seed RNG seed for reproducible intervals
//' 
//' @examples
//' genome <- tibble::frame_data(
//'   ~chrom,  ~size,
//'   "chr1",  10000000,
//'   "chr2",  50000000,
//'   "chr3",  60000000,
//'   "chrX",  5000000
//' ) 
//' 
//' bed_random_cpp(genome, length = 100, n = 1000)
//' 
// [[Rcpp::export]]
DataFrame bed_random_cpp(DataFrame genome, int length, int n, unsigned int seed = 0) {
 
  std::vector<std::string> chroms = as<std::string>(genome["chrom"]) ;
  std::vector<int> sizes = as<int>(genome["size"]) ;
  
  typedef boost::mt19937 RNGtype ;
  RNGtype rng(seed) ;
  
  // chrom RNG 
  int nchrom = chroms.size() ;
  
  boost::uniform_int<> chrom_dist(1, nchrom) ;
  boost::variate_generator< RNGtype, boost::uniform_int<> > chrom_rng(rng, chrom_dist) ;
  
  // size RNGs
  std::vector< boost::variate_generator< RNGtype, boost::uniform_int <> > > size_rngs ;
  
  for (int i=0; i<nchrom; ++i) {
    int size = sizes[i] ;
    // sub length to avoid off-chrom coordinates
    boost::uniform_int<> size_dist(1, size - length) ;
    boost::variate_generator< RNGtype, boost::uniform_int<> > size_rng(rng, size_dist) ;
    size_rngs.push_back(size_rng) ;
  }
 
  std::vector<std::string> rand_chroms ;
  std::vector<int> rand_starts ;
  
  for (int i=0; i<n; ++i) {
  //   
     int chrom_idx = chrom_rng() ;
  //   
     boost::variate_generator< RNGtype,
                               boost::uniform_int<> > size_rng = size_rngs[chrom_idx] ;
  //   
     rand_chroms.push_back(chroms[chrom_idx]) ;
  //   rand_starts.push_back( size_rng() ) ;
  }
  
  std::vector<int> rand_ends ;
  
  return DataFrame::create( Named("chrom") = rand_chroms,
                            Named("start") = rand_starts,
                            Named("end") = rand_ends) ;
  
}
