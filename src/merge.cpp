#include <Rcpp.h>
using namespace Rcpp ;

//[[Rcpp::depends(dplyr,BH)]]
#include <dplyr.h>
using namespace dplyr ;

std::string interval_id(std::string const& chrom, int const& start, int const& end) {
  std::ostringstream id ;
  id << chrom << ':' << start << '-' << end ;
  return id.str() ;
}
int interval_overlap(int const& start_x, int const& end_x, int const& start_y, int const& end_y) {
  return(std::min(end_x, end_y) - std::max(start_x, start_y)) ;    
}

//[[Rcpp::export]]
DataFrame merge_impl(GroupedDataFrame gdf, int max_dist = 0) {
  
  int ng = gdf.ngroups() ;
  
  DataFrame df = gdf.data() ;
  
  int nr = df.nrows() ;
  int nc = df.size() ;
  
  CharacterVector ids(nr) ;
  IntegerVector overlaps(nr) ;
  
  CharacterVector chroms = df["chrom"] ; 
  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;
 
  GroupedDataFrame::group_iterator git = gdf.group_begin() ;
  for(int i=0; i<ng; i++, ++git) {
    
    SlicingIndex indices = *git ;
    int ni = indices.size() ;
   
    int last_start = 0 ;
    int last_end = 0 ;
    std::string last_id ;
  
    for(int j=0; j<ni; j++) {
      int idx = indices[j] ;
     
      std::string chrom = as<std::string>(chroms[idx]) ; 
      int start = starts[idx] ;
      int end = ends[idx] ;
     
      std::string id = interval_id(chrom, start, end) ; 
     
      int overlap = interval_overlap(start, end, last_start, last_end) ;
      overlaps[idx] = overlap ;
     
      if ( overlap > 0 ) {
        ids[idx] = last_id ;
      } else if ( overlap <= 0 && std::abs(overlap) <= max_dist ) {
        ids[idx] = last_id ;
      } else {
        ids[idx] = id ;
        last_id = id ;
      }
      
      last_end = end ; last_start = start ;
    }
  }
  
  // add two new columns, hashes and overlaps
  List out(nc + 2) ;
  CharacterVector onames = df.attr("names") ;
  CharacterVector names( nc + 2 ) ;
  
  for( int i=0; i<nc; i++) {
    out[i] = df[i] ;
    names[i] = onames[i] ;
  }
 
  // add ids 
  out[nc] = ids;
  std::string id_name = ".merge_id" ;
  names[nc] = id_name ;
  
  // add overlaps
  out[nc+1] = overlaps ;
  std::string overlaps_name = ".overlap" ;
  names[nc+1] = overlaps_name;
 
  out.attr("class") = df.attr("class") ;
  out.attr("row.names") = df.attr("row.names") ;
  out.attr("names") = names ;
  
  return out ;
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
//   merge_impl(bed_tbl)
// */
