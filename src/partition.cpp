// partition.cpp
//
// Copyright (C) 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

class IntervalCache {
  // store vector of intervals to partition
  // and max stop of the stored intervals
public:
  ivl_vector_t ivls ;
  int max_stop = 0;
  void clear() {
    ivls.clear();
    max_stop = 0;
  }
  void push_back(ivl_t ivl) {
    ivls.push_back(ivl);
  }
  size_t size() {
    return ivls.size();
  }
};

void partitionIntervals(const IntervalCache& ivl_cache,
                        ivl_vector_t& ivl_result) {
  // convert interval to points
  std::vector<int> ivl_points ;
  for (auto it:ivl_cache.ivls) {
    ivl_points.push_back(it.start) ;
    ivl_points.push_back(it.stop) ;
  }

  // get unique positions
  std::sort(ivl_points.begin(), ivl_points.end()) ;
  ivl_points.erase(std::unique(ivl_points.begin(), ivl_points.end()), ivl_points.end());

  // generate partitions as intervals between
  size_t n_pts = ivl_points.size() ;
  for (size_t i = 0; i < n_pts - 1; i++) {
    ivl_t new_ivl(ivl_points[i], ivl_points[i + 1], ivl_cache.ivls[0].value) ;
    ivl_result.push_back(new_ivl) ;
  }
}

//[[Rcpp::export]]
DataFrame partition_impl(const ValrGroupedDataFrame& gdf,
                         int max_dist = -1) {

  auto ng = gdf.ngroups() ;
  DataFrame df = gdf.data() ;

  ListView idx(gdf.indices()) ;

  std::vector<ivl_t> out_ivls ;
  IntervalCache ivl_cache ;
  for (int i = 0; i < ng; i++) {

    IntegerVector indices ;
    indices = idx[i];
    ivl_vector_t intervals = makeIntervalVector(df, indices);

    // set first interval
    ivl_cache.push_back(intervals[0]) ;
    ivl_cache.max_stop = intervals[0].stop ;
    intervals.erase(intervals.begin()) ;
    for (auto it : intervals) {

      auto max_stop = ivl_cache.max_stop ;
      if (max_stop + max_dist < it.start) {
        // doesn't overlap
        // clear out cache
        if (ivl_cache.ivls.size() > 1) {
          partitionIntervals(ivl_cache,
                             out_ivls) ;
          ivl_cache.clear() ;
        } else {
          out_ivls.push_back(ivl_cache.ivls[0]) ;
          ivl_cache.clear() ;
        }

        ivl_cache.push_back(it) ;
        ivl_cache.max_stop = it.stop ;
      }

      else if (max_stop + max_dist < it.stop) {
        // overlaps
        // update max_stop position and store ivl
        ivl_cache.max_stop = it.stop ;
        ivl_cache.push_back(it) ;
      }

      else {
        // contained interval
        // store interval
        ivl_cache.push_back(it) ;
      }
    }
    // write out final stored interval set
    if (ivl_cache.size() > 1) {
      partitionIntervals(ivl_cache,
                         out_ivls) ;
      ivl_cache.clear() ;
    } else {
      out_ivls.push_back(ivl_cache.ivls[0]) ;
      ivl_cache.clear() ;
    }
  }

  std::vector<int> indices_x ;
  std::vector<int> group_starts ;
  std::vector<int> group_ends ;

  // iterate through vector of partitioned intervals and write to dataframe
  for (auto it : out_ivls) {
    indices_x.push_back(it.value) ;
    group_starts.push_back(it.start) ;
    group_ends.push_back(it.stop) ;
  }

  // subset original dataframe, note that the grouped column values will be correct
  // but non-grouped column values are no longer matched to each interval
  // and are dropped on the R side
  DataFrame subset_x = subset_dataframe(df, indices_x) ;

  subset_x["start"] = group_starts ;
  subset_x["end"] = group_ends ;
  return subset_x ;
}

/*** R
library(dplyr)
bed_tbl <- trbl_interval(
   ~chrom, ~start, ~end, ~value, ~strand,
  'chr1', 100,    500,  10, "+",
  'chr1', 200,    400,  20, "-",
  'chr1', 300,    550,  30, "+",
  'chr1', 550,    575,   2, "+",
  'chr1', 800,    900,   5, "+" ) %>%
  group_by(chrom)

partition_impl(bed_tbl) %>% as_data_frame
*/
