// closest.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void closest_grouped(ivl_vector_t& vx, ivl_vector_t& vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes, std::vector<int>& distance_sizes) {

  ivl_tree_t tree_y(vy) ;

  std::pair<int, ivl_vector_t> min_dist;
  // initiatialize maximum left and right distances to minimize for closest
  int max_end = std::max(vx.back().stop, vy.back().stop) ;

  for (auto const& vx_it : vx) {
    ivl_vector_t closest ;
    ivl_vector_t closest_ivls ;

    min_dist = std::make_pair(max_end, closest_ivls) ;
    tree_y.findClosest(vx_it.start, vx_it.stop, closest, min_dist) ;

    for (auto const& ov_it : closest) {

      auto overlap = intervalOverlap(vx_it, ov_it) ;

      if (overlap > 0) {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
        overlap_sizes.push_back(overlap < 0 ? -overlap : overlap) ;
        distance_sizes.push_back(0);
      } else if (ov_it.start > vx_it.stop) {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
        overlap_sizes.push_back(0) ;
        distance_sizes.push_back(-overlap);
      } else {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
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

  DataFrame df_x = x.data() ;
  DataFrame df_y = y.data() ;

  // for subsetting / return df
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  std::vector<int> overlap_sizes ;
  std::vector<int> distance_sizes ;

  // preallocate vector to save some time
  indices_x.reserve(2e6) ;
  indices_y.reserve(2e6) ;
  overlap_sizes.reserve(2e6) ;
  distance_sizes.reserve(2e6) ;

  // set up interval trees for each chromosome and apply closest_grouped
  GroupApply(x, y, closest_grouped, std::ref(indices_x), std::ref(indices_y),
             std::ref(overlap_sizes), std::ref(distance_sizes));

  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(df_y, names(df_y)).subset(indices_y, "data.frame");

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, suffix_x, false) ;

  // y names, data
  out.add_df(subset_y, suffix_y, true) ;

  // overlaps and distances
  out.add_vec(".overlap", wrap(overlap_sizes)) ;
  out.add_vec(".dist", wrap(distance_sizes)) ;

  auto nrows = subset_x.nrows() ;
  auto res = out.format_df(nrows) ;
  return res ;

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
