#include "valr.h"

void reldist_grouped(intervalVector& vx, intervalVector& vy,
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
  
  // iterate through x intervals and calculate reldist using a binary search
  for (vx_it = vx.begin(); vx_it != vx.end(); ++vx_it ){
    int midpoint = (vx_it->start + vx_it->stop) / 2 ; 
    low_it = std::lower_bound(ref_midpoints.begin(), ref_midpoints.end(), midpoint) ;
    
    low_idx = low_it - ref_midpoints.begin() ;
    
    // drop intervals at start and end which have no reldist
    if (low_idx == 0 ||
        low_idx == ref_midpoints.size()) {
      continue ; 
    }

    // get index below and above  
    low_idx = low_idx - 1;
    upper_idx = low_idx + 1 ;
    
    int left = ref_midpoints[low_idx] ; 
    int right = ref_midpoints[upper_idx] ;
    
    std::size_t dist_l = abs(midpoint - left) ;
    std::size_t dist_r = abs(midpoint - right) ;
    
    //calc relative distance
    float reldist = (float) std::min(dist_l, dist_r) / float(right - left) ;
    
    rel_distances.push_back(reldist) ; 
    indices_x.push_back(vx_it->value) ;
      
  }
}

//[[Rcpp::export]]
DataFrame reldist_impl(GroupedDataFrame x, GroupedDataFrame y) {
  
  std::vector<float> rel_distances ; 
  std::vector<int> indices_x ;
  
  DataFrame df_x = x.data() ;
  DataFrame df_y = y.data() ;
  
  int ng_x = x.ngroups() ;
  int ng_y = y.ngroups() ;
    
  // get labels info for grouping
  DataFrame labels_x(df_x.attr("labels")); 
  DataFrame labels_y(df_y.attr("labels")); 
    
  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for( int nx=0; nx<ng_x; nx++, ++git_x){
      
    SlicingIndex gi_x = *git_x ;
    
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for( int ny=0; ny<ng_y; ny++, ++git_y) {
        
      SlicingIndex gi_y = *git_y ;
        
        // make sure that x and y groups are the same
      bool same_groups = compareDataFrameRows(labels_x, labels_y, nx, ny); 
        
      if(same_groups){
          
        intervalVector vx = makeIntervalVector(df_x, gi_x) ;
        intervalVector vy = makeIntervalVector(df_y, gi_y) ;
          
        reldist_grouped(vx, vy, indices_x, rel_distances);
        }
      } 
    }

  
  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_x, "data.frame");
  
  int ncol_x = subset_x.size() ;
  
  CharacterVector names(ncol_x + 1) ;
  CharacterVector names_x = subset_x.attr("names") ;
  
  List out(ncol_x + 1) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    names[i] = names_x[i] ;
    out[i] = subset_x[i] ;
  }
  out[ncol_x] = rel_distances ;
  names[ncol_x] = "reldist" ;
  
  out.attr("names") = names ; 
  out.attr("class") = classes_not_grouped() ;
  int nrows = subset_x.nrows() ; 
  set_rownames(out, nrows) ;
  
  return out ; 
  
}

/***R
library(dplyr)
  x <- tibble::tribble(
      ~chrom, ~start, ~end,
      "chr1", 5,    15,
      "chr1", 50, 150, 
      "chr2", 1000,   2000,
      "chr3", 3000, 4000
  ) %>% group_by(chrom)
  
  y <- tibble::tribble(
      ~chrom, ~start, ~end,
      "chr1", 25,    125,
      "chr1", 150,    250,
      "chr1", 550,    580,
      "chr2", 1,   1000,
      "chr2", 2000, 3000
  ) %>% group_by(chrom)
  


reldist_impl(x, y)
  
  */
