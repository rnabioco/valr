#include "Rbedtools.h"

void intersect_group(intervalVector vx, intervalVector vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes) {
  
  intervalTree tree_y(vy) ;
  intervalVector overlaps ;

  intervalVector::const_iterator it ; 
  for(it = vx.begin(); it != vx.end(); ++it) {
    
    tree_y.findOverlapping(it->start, it->stop, overlaps) ;
    
    // store current intervals 
    intervalVector::const_iterator oit ; 
    for(oit = overlaps.begin(); oit != overlaps.end(); ++oit) {
      
      int overlap_size = intervalOverlap(*it, *oit) ;
      overlap_sizes.push_back(overlap_size) ;
      
      indices_x.push_back(it->value) ;
      indices_y.push_back(oit->value) ;
    }
    
    overlaps.clear() ;
  } 
}


//' @rdname bed_intersect
//' 
//[[Rcpp::export]]
DataFrame intersect_impl(GroupedDataFrame x, GroupedDataFrame y,
                         std::string suffix_x = ".x", std::string suffix_y = ".y") {
  
  // indices for subsetting 
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  // overlap sizes
  std::vector<int> overlap_sizes ;
  
  auto data_x = x.data() ;
  auto data_y = y.data() ;
  
  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;
  
  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for( int nx=0; nx<ng_x; nx++, ++git_x){
    
    SlicingIndex gi_x = *git_x ;
    auto group_x = gi_x.group() ;
    
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for( int ny=0; ny<ng_y; ny++, ++git_y) {
      
      SlicingIndex gi_y = *git_y ;
      auto group_y = gi_y.group() ;
      
      if( group_x == group_y ) {
        
        intervalVector vx = makeIntervalVector(data_x, gi_x) ;
        intervalVector vy = makeIntervalVector(data_y, gi_y) ;
       
        intersect_group(vx, vy, indices_x, indices_y, overlap_sizes) ; 
      }
    } 
  }
  
  DataFrame subset_x = DataFrameSubsetVisitors(data_x, names(data_x)).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(data_y, names(data_y)).subset(indices_y, "data.frame");
  
  auto ncol_x = subset_x.size() ;
  auto ncol_y = subset_y.size() ;
  
  CharacterVector names(ncol_x + ncol_y) ;
  CharacterVector names_x = subset_x.attr("names") ;
  CharacterVector names_y = subset_y.attr("names") ;

  // replacing y chrom with overlap, same number of cols 
  List out(ncol_x + ncol_y) ;
  
  // x names, data
  for( int i=0; i<ncol_x; i++) {
    auto name_x = as<std::string>(names_x[i]) ;
    if (name_x == "start" || name_x == "end") {
      name_x = name_x + suffix_x ;
    } 
    names[i] = name_x ;
    out[i] = subset_x[i] ;
  }
  
  // y names, data
  for( int i=0; i<ncol_y; i++) {
    auto name_y = as<std::string>(names_y[i]) ;
    
    if (name_y == "chrom") continue ;
    
    if (name_y == "start" || name_y == "end") {
      name_y = name_y + suffix_y ;
    } 
    
    names[i+ncol_x-1] = name_y ;
    out[i+ncol_x-1] = subset_y[i] ;
  }
  
  // overlaps 
  out[ncol_x + ncol_y - 1] = overlap_sizes ;
  names[ncol_x + ncol_y - 1] = ".overlap" ;
  
  out.attr("names") = names ; 
  out.attr("class") = classes_not_grouped() ;
  auto nrows = subset_x.nrows() ; 
  set_rownames(out, nrows) ;
  
  return out ; 
  
}  

/***R
library(Rbedtools)
library(dplyr)

genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7
)

n <- 1e5
x <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)
y <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)

library(microbenchmark)
microbenchmark(
  intersect_impl(x, y)
)
*/
