// closest.cpp
//
// Copyright (C) 2016 - 2022 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"


void findClosest(const IntervalTree<int, int>& tree, int start, int stop,
                                  ivl_vector_t& closest,
                                  std::pair<int, ivl_vector_t>& min_dist) {
  findClosestIvls(tree, start, stop, closest, min_dist);
}


// template<class T, typename K>
// void findClosest(const IntervalTree<T, K>& tree, int start, int stop,
//                  ivl_vector_t& closest,
//                  std::pair<int, ivl_vector_t>& min_dist) {
//
//   ivl_vector_t intervals = tree.intervals;
//   typedef ivl_t interval;
//   if (!intervals.empty() && !(stop < intervals.front().start)) {
//     for (typename ivl_vector_t::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
//       const interval& interval = *i;
//       if (interval.stop > start && interval.start < stop) {
//         // adjacent intervals are considered non-overlappping
//         closest.push_back(interval);
//       } else if (stop <= interval.start) {
//         // find distance on left
//         int ivl_dist_l = interval.start - stop ;
//         // if smaller than best min dist found update pair with dist and intervals
//         if (ivl_dist_l < min_dist.first) {
//           min_dist.first = ivl_dist_l ;
//           min_dist.second.clear() ;
//           min_dist.second.push_back(interval) ;
//         } else if (ivl_dist_l == min_dist.first) {
//           // if same dist append intervals
//           min_dist.second.push_back(interval) ;
//         }
//       } else if (start >= interval.stop) {
//         // find distance on right
//         int ivl_dist_r = start - interval.stop ;
//         // if smaller than best min dist found update pair with dist and intervals
//         if (ivl_dist_r < min_dist.first) {
//           min_dist.first = ivl_dist_r ;
//           min_dist.second.clear() ;
//           min_dist.second.push_back(interval) ;
//         } else if (ivl_dist_r == min_dist.first) {
//           // if same dist append interval
//           min_dist.second.push_back(interval) ;
//         }
//       }
//     }
//   }  else if (!intervals.empty()  && (stop <= intervals.front().start)) {
//     for (typename ivl_vector_t::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
//       const interval& interval = *i;
//       if (interval.start > intervals.front().start) {
//         continue ;
//       } else {
//         // find distance on left
//         int ivl_dist_l = interval.start - stop ;
//         // if smaller than best min dist found update pair with dist and intervals
//         if (ivl_dist_l <= min_dist.first) {
//           min_dist.first = ivl_dist_l;
//           min_dist.second.clear() ;
//           min_dist.second.push_back(interval) ;
//         } else if (ivl_dist_l == min_dist.first) {
//           // if same dist append intervals
//           min_dist.second.push_back(interval) ;
//         }
//       }
//     }
//   }
//
//
//   if (tree.left && start <= tree.center) {
//     findClosest(*tree.left, start, stop, closest, min_dist);
//   }
//
//   if (tree.right && stop >= tree.center) {
//     findClosest(*tree.right, start, stop,closest, min_dist);
//   }
//
//   // Finally report all of the non-overlapping closest intervals, only if at a left_node
//   if (!(tree.right && tree.left)) {
//     closest.insert(closest.end(), min_dist.second.begin(), min_dist.second.end())  ;
//   }
// }

void closest_grouped(const ivl_vector_t& vx, ivl_vector_t& vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes, std::vector<int>& distance_sizes) {

  std::pair<int, ivl_vector_t> min_dist;

  // initialize maximum left and right distances to minimize for closest
  int max_end = std::max(vx.back().stop, vy.back().stop) ;

  // vy should not be referenced after std::move
  ivl_tree_t tree_y(std::move(vy)) ;

  for (auto const& vx_it : vx) {
    ivl_vector_t closest ;
    ivl_vector_t closest_ivls ;

    min_dist = std::make_pair(max_end, closest_ivls) ;
    findClosest(tree_y, vx_it.start, vx_it.stop, closest, min_dist) ;

    for (auto const& ov_it : closest) {

      auto overlap = intervalOverlap(vx_it, ov_it) ;

      if (overlap > 0) {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
        overlap_sizes.push_back(overlap < 0 ? -overlap : overlap) ;
        distance_sizes.push_back(0);
      } else if (ov_it.start >= vx_it.stop) {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
        overlap_sizes.push_back(0) ;
        distance_sizes.push_back(-(overlap - 1));
      } else {
        indices_x.push_back(vx_it.value) ;
        indices_y.push_back(ov_it.value) ;
        overlap_sizes.push_back(0) ;
        distance_sizes.push_back(overlap - 1);
      }

    }
    closest.clear() ;
  }
}

//[[Rcpp::export]]
DataFrame closest_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y,
                       IntegerVector grp_idx_x,
                       IntegerVector grp_idx_y,
                       const std::string& suffix_x,
                       const std::string& suffix_y) {

  DataFrame df_x = x.data() ;
  DataFrame df_y = y.data() ;

  // for subsetting / return df
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  std::vector<int> overlap_sizes ;
  std::vector<int> distance_sizes ;

  // set up interval trees for each chromosome and apply closest_grouped
  GroupApply(x, y, grp_idx_x, grp_idx_y,
             closest_grouped, std::ref(indices_x), std::ref(indices_y),
             std::ref(overlap_sizes), std::ref(distance_sizes));

  DataFrame subset_x = subset_dataframe(df_x, indices_x) ;
  DataFrame subset_y = subset_dataframe(df_y, indices_y) ;

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

closest_impl(x, y, suffix_x, suffix_y, environment())

*/
