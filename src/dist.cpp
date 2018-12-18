// dist.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void dist_grouped(ivl_vector_t& vx, ivl_vector_t& vy,
                  std::vector<int>& indices_x,
                  std::vector<double>& distances,
                  std::string dist_fxn) {

  // first build sorted vector of y interval midpoints

  std::vector<int> ref_midpoints ;
  for (auto vy_it : vy) {
    int midpoint = (vy_it.start + vy_it.stop) / 2 ;
    ref_midpoints.push_back(midpoint) ;
  }

  std::sort(ref_midpoints.begin(), ref_midpoints.end()) ;

  std::size_t low_idx, upper_idx;

  // iterate through x intervals and calculate dist using a binary search
  for (auto vx_it : vx) {
    int midpoint = (vx_it.start + vx_it.stop) / 2 ;
    auto low_it = std::lower_bound(ref_midpoints.begin(),
                                   ref_midpoints.end(), midpoint) ;

    low_idx = low_it - ref_midpoints.begin() ;

    if (dist_fxn == "absdist") {
      // set up indexes for closest element, handling edge cases at start and end of x ivl vector
      if (low_idx == 0) {
        // no need to continue return absdist
        int dist = ref_midpoints[low_idx] ;
        int absdist = abs(dist - midpoint) ;

        distances.push_back(absdist) ;
        indices_x.push_back(vx_it.value) ;
        continue ;

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

      int dist_l = abs(midpoint - left) ;
      int dist_r = abs(midpoint - right) ;

      //calc relative distance
      int absdist = std::min(dist_l, dist_r);

      distances.push_back(absdist) ;
      indices_x.push_back(vx_it.value) ;

    } else {  // reldist
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

      int dist_l = abs(midpoint - left) ;
      int dist_r = abs(midpoint - right) ;

      // calc relative distance
      auto reldist = (float) std::min(dist_l, dist_r) / float(right - left) ;

      distances.push_back(reldist) ;
      indices_x.push_back(vx_it.value) ;
    }
  }
}

//[[Rcpp::export]]
DataFrame dist_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y,
                    IntegerVector x_grp_indexes,
                    IntegerVector y_grp_indexes,
                    std::string distcalc) {

  std::vector<double> distances ;
  std::vector<int> indices_x ;

  DataFrame df_x = x.data() ;
  GroupApply(x, y, x_grp_indexes, y_grp_indexes, dist_grouped, std::ref(indices_x), std::ref(distances), std::ref(distcalc));

  DataFrame subset_x = subset_dataframe(df_x, indices_x) ;

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, false) ;

  // distances
  std::string distname = distcalc == "absdist" ? ".absdist" : ".reldist" ;
  out.add_vec(distname, wrap(distances)) ;

  auto nrows = subset_x.nrows() ;
  auto res = out.format_df(nrows) ;
  return res ;

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
  devtools::load_all()
  dist_impl(x, y, distcalc = "absdist")
  dist_impl(x, y, distcalc = "reldist")
  */
