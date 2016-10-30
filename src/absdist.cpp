#include "valr.h"

void absdist_grouped(intervalVector& vx, intervalVector& vy,
                     std::vector<int>& indices_x,
                     std::vector<float>& rel_distances){
  
  // first build sorted vector of y interval midpoints 
  
  std::vector<int> ref_midpoints ; 
  intervalVector::const_iterator vy_it ;
  for (vy_it = vy.begin(); vy_it != vy.end(); ++vy_it ){
    int midpoint = (vy_it->start + vy_it->stop) / 2 ; 
    ref_midpoints.push_back(midpoint) ;
  }
  
  std::sort(ref_midpoints.begin(), ref_midpoints.end()) ;
  
  std::vector<int>::iterator low_it;
  std::size_t low_idx, upper_idx;
  intervalVector::const_iterator vx_it ;
  
  // iterate through x intervals and calculate absdist using a binary search
  for (vx_it = vx.begin(); vx_it != vx.end(); ++vx_it ){
    int midpoint = (vx_it->start + vx_it->stop) / 2 ; 
    low_it = std::lower_bound(ref_midpoints.begin(), ref_midpoints.end(), midpoint) ;
    
    low_idx = low_it - ref_midpoints.begin() ;
    
    // set up indexes for closest element, handling edge cases at start and end of x ivl vector
    if (low_idx == 0) {
      // just search for closest higher ivl, but do it twice to prevent indexing errors 
      upper_idx = low_idx + 1 ; 
    } else if (low_idx == ref_midpoints.size()) {
      // just search for closest lower ivl
      low_idx = low_idx - 1;
      upper_idx = low_idx;
    } else {
      // search either
      // get index below and above  
      low_idx = low_idx - 1;
      upper_idx = low_idx + 1 ;
    }
    
    int left = ref_midpoints[low_idx] ; 
    int right = ref_midpoints[upper_idx] ;
    
    std::size_t dist_l = abs(midpoint - left) ;
    std::size_t dist_r = abs(midpoint - right) ;
    
    //calc relative distance
    std::size_t absdist = std::min(dist_l, dist_r);
    
    rel_distances.push_back(absdist) ; 
    indices_x.push_back(vx_it->value) ;
    
  }
}

//[[Rcpp::export]]
DataFrame absdist_impl(GroupedDataFrame x, GroupedDataFrame y) {
  
  std::vector<float> rel_distances ; 
  std::vector<int> indices_x ;
  
  DataFrame df_x = x.data() ;
  PairedGroupApply(x, y, absdist_grouped, std::ref(indices_x), std::ref(rel_distances)); 
  
  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_x, "data.frame");
  
  auto ncol_x = subset_x.size() ;
  
  CharacterVector names(ncol_x + 1) ;
  CharacterVector names_x = subset_x.attr("names") ;
  
  List out(ncol_x + 1) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    names[i] = names_x[i] ;
    out[i] = subset_x[i] ;
  }
  out[ncol_x] = rel_distances ;
  names[ncol_x] = ".absdist" ;
  
  out.attr("names") = names ; 
  out.attr("class") = classes_not_grouped() ;
  auto nrows = subset_x.nrows() ; 
  set_rownames(out, nrows) ;
  
  return out ; 
  
}

/***R
library(dplyr)
x <- tibble::frame_data(
  ~chrom, ~start, ~end,
  "chr1", 5,    15,
  "chr1", 50, 150, 
  "chr2", 1000,   2000,
  "chr3", 3000, 4000
) %>% group_by(chrom)

y <- tibble::frame_data(
  ~chrom, ~start, ~end,
  "chr1", 25,    125,
  "chr1", 150,    250,
  "chr1", 550,    580,
  "chr2", 1,   1000,
  "chr2", 2000, 3000
) %>% group_by(chrom)



absdist_impl(x, y)

*/
