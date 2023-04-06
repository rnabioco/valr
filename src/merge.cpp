// merge.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

DataFrame collapseMergedIntervals(const ValrGroupedDataFrame& gdf,
                                  int max_dist = 0) {

  auto ng = gdf.ngroups() ;
  DataFrame df = gdf.data() ;

  ListView idx(gdf.indices()) ;

  // approach from http://www.geeksforgeeks.org/merging-intervals/
  std::vector<ivl_t> s ;
  for (int i = 0; i < ng; i++) {

    IntegerVector indices ;
    indices = idx[i];

    ivl_vector_t intervals = makeIntervalVector(df, indices);

    // set first interval
    s.push_back(intervals[0]) ;
    intervals.erase(intervals.begin()) ;
    for (auto it : intervals) {

      auto top = s.back() ;
      if (top.stop + max_dist < it.start) {
        // no overlap push to stack
        s.push_back(it) ;
      }

      else if (top.stop < it.stop) {
        // overlaps and need to update stack top position
        top.stop = it.stop ;
        s.pop_back() ;
        s.push_back(top) ;
      }
    }
  }

  std::vector<int> indices_x ;
  std::vector<int> group_starts ;
  std::vector<int> group_ends ;
  // interate through vector of merged intervals and write to dataframe
  for (auto it : s) {
    indices_x.push_back(it.value) ;
    group_starts.push_back(it.start) ;
    group_ends.push_back(it.stop) ;
  }

  DataFrame subset_x = subset_dataframe(df, indices_x) ;

  subset_x["start"] = group_starts ;
  subset_x["end"] = group_ends ;
  return subset_x ;
}

DataFrame clusterMergedIntervals(const ValrGroupedDataFrame& gdf,
                                 int max_dist = 0) {

  auto ng = gdf.ngroups() ;
  DataFrame df = gdf.data() ;

  auto nr = df.nrows() ;

  IntegerVector ids(nr) ;      // store ids
  std::size_t cluster_id = 0;  // counter for cluster id

  ListView idx(gdf.indices()) ;

  for (int i = 0; i < ng; i++) {

    IntegerVector indices ;
    indices = idx[i];
    int ni = indices.size();

    ivl_vector_t intervals = makeIntervalVector(df, indices);

    // store an interval to ensure first interval maintained
    ivl_t last_interval = ivl_t(0, 0, 0) ;
    ivl_t top = last_interval;

    for (int j = 0; j < ni ; j++) {
      ivl_t it = intervals[j];
      // get index to store cluster ids in vector
      auto idx = it.value ;
      last_interval = it ; // set interval to compare

      if ( j == 0 || top.stop + max_dist < it.start) {
        // no overlap, update interval and cluster id
        top = it ;
        cluster_id++ ;
        ids[idx] = cluster_id ;
      } else {
        // overlaps or contained in stack top ivl
        if(it.stop > top.stop){
          top.stop = it.stop ; // update end position
        }
        ids[idx] = cluster_id ;
      }
    }
  }

  DataFrameBuilder out;
  // x names, data
  out.add_df(df, false) ;

  // ids
  out.names.push_back(".id_merge");
  out.data.push_back(ids);

  auto nrows = df.nrows() ;
  auto res = out.format_df(nrows) ;
  return res ;

}

//[[Rcpp::export]]
DataFrame merge_impl(ValrGroupedDataFrame gdf,
                     int max_dist = 0,
                     bool collapse = true) {

  if (!collapse) {
    // return a cluster id per input interval
    DataFrame out = clusterMergedIntervals(gdf, max_dist) ;
    return out ;
  } else {
    // return only merged intervals
    DataFrame out = collapseMergedIntervals(gdf, max_dist) ;
    return out ;
  }
}

/*** R
  library(dplyr)
  bed_tbl <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 1,      50,   1,
    "chr1", 100,    200,  2,
    "chr1", 150,    250,  3,
    "chr1", 175,    225,  3.5,
    "chr2", 1,      25,   4,
    "chr2", 200,    400,  5,
    "chr2", 400,    500,  6,
    "chr2", 450,    550,  7,
    "chr3", 450,    550,  8,
    "chr3", 500,    600,  9
  ) %>% group_by(chrom)

  merge_impl(bed_tbl)
*/
