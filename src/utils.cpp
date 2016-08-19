#include "valr.h"

genome_map_t makeChromSizes(DataFrame genome) {
  
  genome_map_t chrom_sizes ; 
  
  CharacterVector chroms = genome["chrom"] ;
  IntegerVector sizes = genome["size"] ;
  
  int nchrom = genome.nrows() ;
  for( int i=0; i<nchrom; ++i) {
    std::string chrom = as<std::string>(chroms[i]) ;
    int size = sizes[i] ;
    chrom_sizes.insert({chrom, size}) ;
  }
 
  return chrom_sizes ; 
}
 
// the value field of intervals in the returned vector correspond to the index
// of the interval in the original dataframe (i.e., the values of the
// SlicingIndex)
intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) {
  
  intervalVector intervals ;
  
  IntegerVector starts = df["start"] ;
  IntegerVector ends   = df["end"] ;
  
  int size = si.size() ;
  
  for( int i=0; i<size; ++i) {
    int j = si[i] ;
    intervals.push_back(interval_t(starts[j], ends[j], j)) ;
  }
  return intervals ;
}

icl_interval_set_t makeIclIntervalSet(DataFrame df, SlicingIndex indices) {

  icl_interval_set_t iset ;
  
  IntegerVector starts = df["start"] ;
  IntegerVector ends   = df["end"] ;
  
  int size = indices.size() ;
  
  for( int i=0; i<size; ++i) {
    int j = indices[i] ;
    icl_interval_t intvl = icl_interval_t::closed(starts[j], ends[j]) ;
    iset.add(intvl) ;
  }

  return iset ;
}

