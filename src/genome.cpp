// genome.cpp
// 
// get genome information

#include <Rcpp.h>
using namespace Rcpp ;

std::map<std::string, int> create_genome(Rcpp::DataFrame genome) {
  
  std::vector <std::string> genome_chrom_v = 
    Rcpp::as< std::vector <std::string> >(genome["chrom"]); 
  std::vector <int> genome_size_v = 
    Rcpp::as< std::vector <int> >(genome["size"]); 
 
  std::map<std::string, int> out ;
  
  for (unsigned i = 0; i < genome_chrom_v.size(); ++i) {
    out[ genome_chrom_v[i] ] = genome_size_v[i] ; 
  }
 
  return out ; 
}
