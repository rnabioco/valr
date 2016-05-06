#include "valr.h"

//' @rdname bed_subtract
//'
//[[Rcpp::export]]
DataFrame subtract_impl(GroupedDataFrame gdf_x, GroupedDataFrame gdf_y) {

  std::vector<std::string> chrom_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
 
  auto ng_x = gdf_x.ngroups() ;
  auto ng_y = gdf_y.ngroups() ;
  
  DataFrame df_x = gdf_x.data() ;
  DataFrame df_y = gdf_y.data() ; 
  
  CharacterVector chroms = df_x["chrom"] ;
  
  GroupedDataFrame::group_iterator git_x = gdf_x.group_begin() ;
  
  for(int nx=0; nx<ng_x; nx++, ++git_x) {
    
    SlicingIndex indices_x = *git_x ; 
    auto group_x = indices_x.group() ;
    
    GroupedDataFrame::group_iterator git_y = gdf_y.group_begin() ;
    for(int ny=0; ny<ng_y; ny++, ++git_y) {
  
      SlicingIndex indices_y = *git_y ; 
      auto group_y = indices_y.group() ;
      
      if ( group_x == group_y ) {

  	    icl_interval_set_t interval_set_x = makeIclIntervalSet(df_x, indices_x) ; 
  	    icl_interval_set_t interval_set_y = makeIclIntervalSet(df_y, indices_y) ;
  	  
  	    // subtract the sets, `-` is overloaded in boost::icl
  	    icl_interval_set_t interval_sub =  interval_set_x - interval_set_y ;
  	    
  	    if (interval_sub.empty()) continue ;
  	   
  	    // get chrom name based on first index in indices_x
  	    // XXX it would be a lot nicer to fectch this from the current data 
  	    // in indices_y via symbol or label but that doesn't seem to work.
  	    std::string chrom = as<std::string>(chroms[indices_x[0]]) ;
  	    
        icl_interval_set_t::iterator it ;
        for( it = interval_sub.begin(); it != interval_sub.end(); ++it) {
          
          chrom_out.push_back(chrom) ;
          starts_out.push_back(it->lower()) ;
          ends_out.push_back(it->upper()) ;
        }
        
      }
    }
	}
  
  return DataFrame::create( Named("chrom") = chrom_out,
                            Named("start") = starts_out,
                            Named("end") = ends_out) ;
}

/***R
library(dplyr)
library(valr)
  
genome <- tibble::frame_data(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7
)
  
n <- 1e4
x <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)
y <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)

subtract_impl(x, y) %>% as_data_frame()

x <- tibble::frame_data(
  ~chrom, ~start, ~end,
  "chr1", 100,    200
) %>% group_by(chrom)

y <- tibble::frame_data(
  ~chrom, ~start, ~end,
  "chr1", 1000,    2000
) %>% group_by(chrom)

subtract_impl(x, y)

*/
