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
        // no overlap push to stack and get new id
        s.push_back(it) ;
      }

      else if (top.stop + max_dist < it.stop) {
        // overlaps and need to update stack top position
        // do not update id
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
  IntegerVector overlaps(nr) ; // store overlap values

  std::size_t cluster_id = 0;  //store counter for cluster id

  ListView idx(gdf.indices()) ;

  // approach from http://www.geeksforgeeks.org/merging-intervals/
  for (int i = 0; i < ng; i++) {

    IntegerVector indices ;
    indices = idx[i];
    int ni = indices.size();

    ivl_vector_t intervals = makeIntervalVector(df, indices);

    // store a first interval to ensure first interval maintained
    ivl_t last_interval = ivl_t(0, 0, 0) ;
    std::vector<ivl_t> s ;
    s.push_back(last_interval) ;

    for (int j = 0; j < ni ; j++) {
      ivl_t it = intervals[j];
      // get index to store cluster ids in vector
      auto idx = it.value ;

      // get overlap distances and assign at proper index
      int overlap = intervalOverlap(it, last_interval) ;
      overlaps[idx] = overlap ;

      last_interval = it ; // set interval to compare
      auto top = s.back() ; // last good interval
      if ( j == 0 || top.stop + max_dist < it.start) {
        // no overlap push to end of vector and get new id
        s.push_back(it) ;
        cluster_id++ ;
        ids[idx] = cluster_id ; // assign cluster id at proper index
      }

      else if (top.stop + max_dist < it.stop) {
        // overlaps and need to update stack top position
        // do not update id
        top.stop = it.stop ; // update end position
        s.pop_back() ; // remove best end
        s.push_back(top) ; // update end interval
        ids[idx] = cluster_id ; // assign cluster id at proper index
      }

      else {
        // overlaps but contained in stack top ivl
        // do not update id
        ids[idx] = cluster_id ;
      }
    }
  }

  DataFrameBuilder out;
  // x names, data
  out.add_df(df, false) ;

  // ids and overlaps
  out.names.push_back(".id_merge");
  out.data.push_back(ids);

  out.names.push_back(".overlap_merge");
  out.data.push_back(overlaps);

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

  merge_impl(bed_tbl) %>% as_data_frame
*/
