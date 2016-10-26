#include "valr.h"

void closest_grouped(intervalVector& vx, intervalVector& vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes, std::vector<int>& distance_sizes) {
  
  intervalTree tree_y(vy) ;
  
  std::pair<int, intervalVector> min_dist;
  // initiatialize maximum left and right distances to minimize for closest
  int max_end = std::max(vx.back().stop, vy.back().stop) ;
  
  intervalVector::const_iterator vx_it ;
  for(vx_it = vx.begin(); vx_it != vx.end(); ++vx_it ) {
    intervalVector closest ;
    intervalVector closest_ivls ; 
    
    min_dist = std::make_pair(max_end, closest_ivls) ;
    tree_y.findClosest(vx_it->start, vx_it->stop, closest, min_dist) ;
    
    intervalVector::const_iterator ov_it ;
    for(ov_it = closest.begin(); ov_it != closest.end(); ++ov_it ) {
     
      int overlap = intervalOverlap(*vx_it, *ov_it) ;
      
      if(overlap > 0) {
        indices_x.push_back(vx_it->value) ;
        indices_y.push_back(ov_it->value) ;
        overlap_sizes.push_back(overlap) ;
        distance_sizes.push_back(0);
      } else if (ov_it->start > vx_it->stop) {
        indices_x.push_back(vx_it->value) ;
        indices_y.push_back(ov_it->value) ;
        overlap_sizes.push_back(0) ;
        distance_sizes.push_back(-overlap);
      } else {
        indices_x.push_back(vx_it->value) ;
        indices_y.push_back(ov_it->value) ;
        overlap_sizes.push_back(0) ;
        distance_sizes.push_back(overlap);
      }

    }
    closest.clear() ; 
  } 
}   

//[[Rcpp::export]]
DataFrame closest_impl(GroupedDataFrame x, GroupedDataFrame y,
                       const std::string& suffix_x, const std::string& suffix_y) {

  // for subsetting / return df
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  std::vector<int> overlap_sizes ;
  std::vector<int> distance_sizes ;
  
  DataFrame data_x = x.data() ;
  DataFrame data_y = y.data() ;
  
  int ng_x = x.ngroups() ;
  int ng_y = y.ngroups() ;
  
  // get labels info for grouping
  DataFrame labels_x(data_x.attr("labels")); 
  DataFrame labels_y(data_y.attr("labels")); 
  
  // set up interval vectors for each group and apply closest_group
  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for( int nx=0; nx<ng_x; nx++, ++git_x){
    
    SlicingIndex gi_x = *git_x ;
    
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for( int ny=0; ny<ng_y; ny++, ++git_y) {
      
      SlicingIndex gi_y = *git_y ;
      
      // make sure that x and y groups are the same
      bool same_groups = compareDataFrameRows(labels_x, labels_y, nx, ny); 
      
      if(same_groups){
        
        intervalVector vx = makeIntervalVector(data_x, gi_x) ;
        intervalVector vy = makeIntervalVector(data_y, gi_y) ;
        
        closest_grouped(vx, vy, indices_x, indices_y,
                        overlap_sizes, distance_sizes); 
      }
    } 
  }
  
  DataFrame subset_x = DataFrameSubsetVisitors(data_x, names(data_x)).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(data_y, names(data_y)).subset(indices_y, "data.frame");
  
  int ncol_x = subset_x.size() ;
  int ncol_y = subset_y.size() ;
  
  CharacterVector names(ncol_x + ncol_y + 1) ;
  CharacterVector names_x = subset_x.attr("names") ;
  CharacterVector names_y = subset_y.attr("names") ;
  
  // replacing y chrom with overlap, and adding distance col
  List out(ncol_x + ncol_y + 1) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    std::string name_x = as<std::string>(names_x[i]) ;
    if (name_x != "chrom") {
      name_x += suffix_x ;
    } 
    names[i] = name_x ;
    out[i] = subset_x[i] ;
  }
  
  // y names, data
  for( int i=0; i<ncol_y; i++) {
    std::string name_y = as<std::string>(names_y[i]) ;
    
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
  int nrows = subset_x.nrows() ; 
  set_rownames(out, nrows) ;
  
  return out ; 
  
}

/***R
library(dplyr)
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 500,    600,
  "chr2", 5000,   6000
) %>% group_by(chrom)

y <- tibble::tribble(
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
