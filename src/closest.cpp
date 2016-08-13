#include "valr.h"

void closest_grouped(intervalVector& vx, intervalVector& vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes, std::vector<int>& distance_sizes) {

  intervalVector::const_iterator vx_it ;
  for(vx_it = vx.begin(); vx_it != vx.end(); ++vx_it ) {
   
    // storage for minimum distance and intervals 
    int left_min, right_min = 0;
    
    interval_t ivl_left  = interval_t(0,0,0);
    interval_t ivl_right = interval_t(0,0,0);
    int left_overlap, right_overlap = 0 ;
    
    intervalVector::const_iterator vy_it ;
    for(vy_it = vy.begin(); vy_it != vy.end(); ++vy_it ) {
     
      auto overlap = intervalOverlap(*vx_it, *vy_it) ;
      
      // store overlapping 
      if(overlap > 0) {
        indices_x.push_back(vx_it->value) ;
        indices_y.push_back(vy_it->value) ;
        overlap_sizes.push_back(overlap) ;
        distance_sizes.push_back(0);
        continue ;
      }
      
      // left y intervals
      if(vy_it->stop <= vx_it->start) {
        if(ivl_left.value == 0 || vy_it->stop > ivl_left.stop) {
          ivl_left = *vy_it ;
          ivl_left.value += 1 ; // increment by 1 to avoid losing first left interval if in SlicingIndex 0
          left_overlap = overlap; 
        }
      }
      
      // right y intervals
      if(vy_it->start >= vx_it->stop) {
                     
        if(ivl_right.value != 0 && vy_it->start > ivl_right.stop) {
          break ;
        } else {
          ivl_right = *vy_it ;
          ivl_right.value += 1 ; // increment by 1 to avoid losing first left interval if in SlicingIndex 0
          right_overlap = overlap; 
        }
      }
      
    } // for y
    
    // store left and right intervals 
    // touching intervals stored as distance +/- 1 (Depending on left/right orientation)
    if (ivl_left.value != 0) {
      ivl_left.value -= 1; // decrease by 1 to return index to correct pos (see line 35)
      indices_x.push_back(vx_it->value) ;
      indices_y.push_back(ivl_left.value) ;
      overlap_sizes.push_back(left_overlap) ;
      distance_sizes.push_back(left_overlap - 1);
    }
    
    if (ivl_right.value != 0) {
      ivl_right.value -= 1; // decrease by 1 to return index to correct pos (see line 47)
      indices_x.push_back(vx_it->value) ; 
      indices_y.push_back(ivl_right.value) ;
      overlap_sizes.push_back(right_overlap) ;
      distance_sizes.push_back(-right_overlap + 1);
    }

  } // for x
}   // void

// XXX the following is verbatim from intersect.cpp except for fxn call  
// XXX and an additional column (distance) is added to output df. should be reused
 
//[[Rcpp::export]]
DataFrame closest_impl(GroupedDataFrame x, GroupedDataFrame y,
                       const std::string& suffix_x, const std::string& suffix_y) {

//  auto ng_x = x.ngroups() ;
//  auto ng_y = y.ngroups() ;

  DataFrame df_x = x.data() ;
  DataFrame df_y = y.data() ;
  
//  auto nr_x = df_x.nrows() ;
//  auto nr_y = df_y.nrows() ;
 
  // for subsetting / return df
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  std::vector<int> overlap_sizes ;
  std::vector<int> distance_sizes ;
  
  // set up interval trees for each chromosome and apply closest_grouped
  chromLoop(x, y, closest_grouped, std::ref(indices_x), std::ref(indices_y), 
            std::ref(overlap_sizes), std::ref(distance_sizes)); 
  
  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(df_y, names(df_y)).subset(indices_y, "data.frame");
  
  auto ncol_x = subset_x.size() ;
  auto ncol_y = subset_y.size() ;
  
  CharacterVector names(ncol_x + ncol_y + 1) ;
  CharacterVector names_x = subset_x.attr("names") ;
  CharacterVector names_y = subset_y.attr("names") ;
  
  // replacing y chrom with overlap, and adding distance col
  List out(ncol_x + ncol_y + 1) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    auto name_x = as<std::string>(names_x[i]) ;
    if (name_x != "chrom") {
      name_x += suffix_x ;
    } 
    names[i] = name_x ;
    out[i] = subset_x[i] ;
  }
  
  // y names, data
  for( int i=0; i<ncol_y; i++) {
    auto name_y = as<std::string>(names_y[i]) ;
    
    if (name_y == "chrom") continue ;
    
    name_y += suffix_y ;
    
    names[i+ncol_x-1] = name_y ;
    out[i+ncol_x-1] = subset_y[i] ;
  }
  
  // overlaps 
  out[ncol_x + ncol_y - 1] = overlap_sizes ;
  names[ncol_x + ncol_y - 1] = ".overlap" ;

  //distances
  out[ncol_x + ncol_y] = distance_sizes ;
  names[ncol_x + ncol_y] = ".distance";
  
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
  "chr1", 500,    600,
  "chr2", 5000,   6000
) %>% group_by(chrom)

y <- tibble::frame_data(
  ~chrom, ~start, ~end,
  "chr1", 100,    200,
  "chr1", 150,    200,
  "chr1", 550,    580,
  "chr2", 7000,   8500
) %>% group_by(chrom)

suffix_x <- '.x'
suffix_y <- '.y'

closest_impl(x, y, suffix_x, suffix_y)
             
*/
