#include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::depends(dplyr,BH)]]
#include <dplyr.h>
using namespace dplyr ;

std::size_t label_hasher(int const& start, int const& end) {
  
  std::size_t seed = 0;
  boost::hash_combine(seed, start) ;
  boost::hash_combine(seed, end) ;
  
  return seed ;
}

int interval_overlap2(int const& start_x, int const& end_x, int const& start_y, int const& end_y) {
  return(std::min(end_x, end_y) - std::max(start_x, start_y)) ;    
}

// [[Rcpp::export]]
CharacterVector merge_labels(GroupedDataFrame gdf, int max_dist = 0) {
  
  int ng = gdf.ngroups() ;
  int nr = gdf.nrows() ;
  
  std::vector<std::size_t> res(nr) ;
  
  DataFrame df = gdf.data() ;
  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;
  
  GroupedDataFrame::group_iterator git = gdf.group_begin() ;
  for(int i=0; i<ng; i++, ++git) {
    
    SlicingIndex indices = *git ;
    int ni = indices.size() ;
   
    int last_start = 0 ;
    int last_end = 0 ;
    std::size_t last_hash = 0 ;
  
    for(int j=0; j<ni; j++) {
      int idx = indices[j] ;
      
      int start = starts[idx] ;
      int end = ends[idx] ;
     
      std::size_t hash = label_hasher(start, end) ;

      if (interval_overlap2(start, end, last_start, last_end) > max_dist) {
        
        res[idx] = last_hash ;
        last_end = end ; last_start = start ;
        continue ;
        
      } else {
        res[idx] = hash ;
      }
      
      last_end = end ; last_start = start ; last_hash = hash ;
    }
  }
  
  return Rcpp::wrap(res) ;
  
}

// /*** R
//   library(dplyr)
//   bed_tbl <- tibble::frame_data(
//     ~chrom, ~start, ~end, ~value,
//     "chr1", 1,      50,   1,
//     "chr1", 100,    200,  2,
//     "chr1", 150,    250,  3,
//     "chr1", 175,    225,  3.5,
//     "chr2", 1,      25,   4,
//     "chr2", 200,    400,  5,
//     "chr2", 400,    500,  6,
//     "chr2", 450,    550,  7,
//     "chr3", 450,    550,  8,
//     "chr3", 500,    600,  9
//   ) %>% group_by(chrom)
//   
//   merge_labels(bed_tbl)
// */
