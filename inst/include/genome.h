#include "valr.h"

typedef std::unordered_map<std::string, int> genome_map_t ;
extern genome_map_t makeChromSizes(DataFrame genome) ;

genome_map_t makeChromSizes(DataFrame genome) {

  genome_map_t chrom_sizes ;

  CharacterVector chroms = genome["chrom"] ;
  IntegerVector sizes = genome["size"] ;

  int nchrom = genome.nrows() ;
  for (int i = 0; i < nchrom; ++i) {
    std::string chrom = as<std::string>(chroms[i]) ;
    int size = sizes[i] ;
    chrom_sizes.insert({chrom, size}) ;
  }

  return chrom_sizes ;
}

