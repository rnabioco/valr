#include "valr.h"

//[[Rcpp::export]]
DataFrame subtract_impl(GroupedDataFrame gdf_x, GroupedDataFrame gdf_y) {

  std::vector<std::string> chrom_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
 
  auto ng_x = gdf_x.ngroups() ;
  auto ng_y = gdf_y.ngroups() ;
  
  DataFrame df_x = gdf_x.data() ;
  DataFrame df_y = gdf_y.data() ; 
  
  CharacterVector chroms_x = df_x["chrom"] ;
  CharacterVector chroms_y = df_y["chrom"] ;
  
  GroupedDataFrame::group_iterator git_x = gdf_x.group_begin() ;
  
  for(int nx=0; nx<ng_x; nx++, ++git_x) {
    
    SlicingIndex indices_x = *git_x ; 
    std::string chrom_x = as<std::string>(chroms_x[indices_x[0]]);
    // keep track of if x chrom is present in y
    bool chrom_x_seen(false);
    
    GroupedDataFrame::group_iterator git_y = gdf_y.group_begin() ;
    for(int ny=0; ny<ng_y; ny++, ++git_y) {
  
      SlicingIndex indices_y = *git_y ; 
      std::string chrom_y = as<std::string>(chroms_y[indices_y[0]]);
      
      if( chrom_x == chrom_y ) {
        chrom_x_seen = true ;
  	    icl_interval_set_t interval_set_x = makeIclIntervalSet(df_x, indices_x) ; 
  	    icl_interval_set_t interval_set_y = makeIclIntervalSet(df_y, indices_y) ;
  	  
  	    // subtract the sets, `-` is overloaded in boost::icl
  	    icl_interval_set_t interval_sub =  interval_set_x - interval_set_y ;
  	    
  	    if (interval_sub.empty()) continue ;
  	    
        icl_interval_set_t::iterator it ;
        for( it = interval_sub.begin(); it != interval_sub.end(); ++it) {
          
          chrom_out.push_back(chrom_x) ;
          starts_out.push_back(it->lower()) ;
          ends_out.push_back(it->upper()) ;
        }
        
      }
    }
    // return x intervals if x chromosome not found in y
    if (chrom_x_seen) {
      continue;
      }
    else {
      DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_x, "data.frame");
      std::vector<std::string> x_chr = subset_x["chrom"] ;
      std::vector<int> x_str = subset_x["start"] ;
      std::vector<int> x_end = subset_x["end"] ;
      
      chrom_out.insert(chrom_out.end(), x_chr.begin(), x_chr.end());
      starts_out.insert(starts_out.end(), x_str.begin(), x_str.end());
      ends_out.insert(ends_out.end(), x_end.begin(), x_end.end());
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
