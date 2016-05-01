#include <dplyr.h>
using namespace dplyr ;
// [[Rcpp::depends(dplyr)]]

#include "Rbedtools.h"

// int interval_overlap(int const& start_x, int const& end_x, int const& start_y, int const& end_y) {
//   return(std::min(end_x, end_y) - std::max(start_x, start_y)) ;    
// }

IntegerVector scan_cache(IntegerVector cache, int idx_x, 
                         IntegerVector& starts_x, IntegerVector& ends_x,
                         IntegerVector& starts_y, IntegerVector& ends_y,
                         IntegerVector& indices_x, IntegerVector& indices_y,
                         IntegerVector& overlaps, int max_dist) {
  
  int start_x = starts_x[idx_x] ;
  int end_x = ends_x[idx_x] ;
  
  IntegerVector temp_cache ;
  
  int nc = cache.size() ;

  for( int i=0; i<nc; i++) {
    
    int idx_y = cache[i] ;
    int start_y = starts_y[idx_y] ;
    int end_y = ends_y[idx_y] ;

    if( start_x <= end_y) {
      temp_cache.push_back(idx_y) ;
  
      int overlap = interval_overlap(start_x, end_x, start_y, end_y) ;
      if( overlap > max_dist ) { 
        indices_x.push_back(idx_x) ;
        indices_y.push_back(idx_y) ;
        overlaps.push_back(overlap) ;
      }
    }        
  }

  return temp_cache ;
}


void intersect_grouped(DataFrame x_data, DataFrame y_data,
                       SlicingIndex gidx_x, SlicingIndex gidx_y,
                       IntegerVector& indices_x, IntegerVector& indices_y,
                       IntegerVector& overlaps, int max_dist) {
 
  IntegerVector cache ;
  
  int size_x = gidx_x.size() ;
  int size_y = gidx_y.size() ;
 
  IntegerVector starts_x = x_data["start"] ;  
  IntegerVector ends_x = x_data["end"] ;  
  IntegerVector starts_y = y_data["start"] ;  
  IntegerVector ends_y = y_data["end"] ;  
 
  for( int i=0; i<size_x; ++i ) {
   
    int idx_x = gidx_x[i] ; 
    
    cache = scan_cache(cache, idx_x,
                       starts_x, ends_x, starts_y, ends_y,
                       indices_x, indices_y,
                       overlaps, max_dist) ;
    
    int start_x = starts_x[idx_x] ;
    int end_x = ends_x[idx_x] ;
   
    for (int j=0; j<size_y; ++j) {
      
      int idx_y = gidx_y[j] ; 
      int start_y = starts_y[idx_y] ;
      int end_y = ends_y[idx_y] ;
      
      if( start_x <= end_y) {
        cache.push_back(idx_y) ;
       
        int overlap = interval_overlap(start_x, end_x, start_y, end_y) ;
        if (overlap > max_dist ) {
          indices_x.push_back(idx_x) ; 
          indices_y.push_back(idx_y) ; 
          overlaps.push_back(overlap) ;
        }
      }
    } 
  }
}

//' @rdname bed_intersect
//' 
// [[Rcpp::export]]
DataFrame intersect_impl(GroupedDataFrame x, GroupedDataFrame y, int max_dist = 0,
                         std::string suffix_x = ".x", std::string suffix_y = ".y") {

  // indices for subsetting 
  IntegerVector indices_x ;
  IntegerVector indices_y ;
  
  IntegerVector overlaps ;
  
  DataFrame data_x = x.data() ;
  DataFrame data_y = y.data() ;
  
  int ng_x = x.ngroups() ;
  int ng_y = y.ngroups() ;
  
  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for( int nx=0; nx<ng_x; nx++, ++git_x){
    
    SlicingIndex gi_x = *git_x ;
    int index_x = gi_x.group() ;
    
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for( int ny=0; ny<ng_y; ny++, ++git_y) {
     
      SlicingIndex gi_y = *git_y ;
      int index_y = gi_y.group() ;

      if( index_x == index_y ) {
        intersect_grouped(data_x, data_y, gi_x, gi_y, indices_x, indices_y, overlaps, max_dist) ; 
      }
    } 
  }
  
  DataFrame subset_x = DataFrameSubsetVisitors(data_x, names(data_x)).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(data_y, names(data_y)).subset(indices_y, "data.frame");
  
  int ncol_x = subset_x.size() ;
  int ncol_y = subset_y.size() ;
  
  CharacterVector names(ncol_x + ncol_y) ;
  CharacterVector names_x = subset_x.attr("names") ;
  CharacterVector names_y = subset_y.attr("names") ;

  // replacing y chrom with overlap, same number of cols
  List out(ncol_x + ncol_y) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    std::string name_x = as<std::string>(names_x[i]) ;
    if (name_x == "chrom") {
      names[i] = name_x ;
    } else {
      names[i] = name_x + suffix_x ;
    }
    out[i] = subset_x[i] ;
  }
  
  // y names, data
  for( int i=0; i<ncol_y; i++) {
    std::string name_y = as<std::string>(names_y[i]) ;
    if (name_y == "chrom") {
      continue ;
    } else {
      names[i+ncol_x-1] = name_y + suffix_y ;
    }
    out[i+ncol_x-1] = subset_y[i] ;
  }
 
  // overlaps 
  out[ncol_x + ncol_y - 1] = overlaps ;
  names[ncol_x + ncol_y - 1] = ".overlap" ;
 
  out.attr("names") = names ; 
  out.attr("class") = classes_not_grouped() ;
  int nrows = subset_x.nrows() ; 
  set_rownames(out, nrows) ;
  
  return out ; 
}

/***R
library(dplyr)
x <- tibble::frame_data(
~chrom, ~start, ~end,
"chr1", 100,    500,
"chr2", 200,    400,
"chr2", 300,    500,
"chr2", 800,    900
) %>% group_by(chrom)

y <- tibble::frame_data(
~chrom, ~start, ~end,
"chr1", 150,    400,
"chr2", 230,    430,
"chr2", 350,    430
) %>% group_by(chrom)
  
intersect_impl(x, y, max_dist = 0)
*/

