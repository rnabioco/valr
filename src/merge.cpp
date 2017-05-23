// merge.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

//[[Rcpp::export]]
DataFrame merge_impl(GroupedDataFrame gdf, int max_dist = 0, bool dots = false ) {

  auto ng = gdf.ngroups() ;

  DataFrame df = gdf.data() ;

  auto nr = df.nrows() ;
  auto nc = df.size() ;

  IntegerVector ids(nr) ;      // store ids
  IntegerVector overlaps(nr) ; // store overlap values

  std::size_t cluster_id = 0;  //store counter for cluster id

  GroupedDataFrame::group_iterator git = gdf.group_begin() ;

  // approach from http://www.geeksforgeeks.org/merging-intervals/

  std::stack<ivl_t> s ;

  for (int i = 0; i < ng; i++, ++git) {

    SlicingIndex indices = *git ;

    ivl_vector_t intervals = makeIntervalVector(df, indices);

    ivl_t last_interval = ivl_t(0, 0, 0) ;

    s.push(intervals[0]) ;
    intervals.erase(intervals.begin()) ;
    for (auto it : intervals) {

      auto idx = it.value ;

      int overlap = intervalOverlap(it, last_interval) ;
      overlaps[idx] = overlap ;
      last_interval = it ;

      auto top = s.top() ;
      if (top.stop + max_dist < it.start) {
        // no overlap push to stack and get new id
        s.push(it) ;
        cluster_id++ ;
        ids[idx] = cluster_id ;
      }

      else if (top.stop + max_dist < it.stop) {
        // overlaps and need to update stack top position
        // do not update id
        top.stop = it.stop ;
        s.pop() ;
        s.push(top) ;
        ids[idx] = cluster_id ;
      }

      else {
        // overlaps but contained in stack top ivl
        // do not update id
        ids[idx] = cluster_id ;
      }
    }
  }

  // store indices to keep if internal merge
  if(!dots){
    std::vector<int> indices_x ;
    std::vector<int> group_starts ;
    std::vector<int> group_ends ;
    while (!s.empty())
    {
      ivl_t t = s.top();
      indices_x.push_back(t.value) ;
      group_starts.push_back(t.start) ;
      group_ends.push_back(t.stop) ;
      s.pop();
    }
    DataFrame subset_x = DataFrameSubsetVisitors(df, names(df)).subset(indices_x, "data.frame");
    subset_x["start"] = group_starts ;
    subset_x["end"] = group_ends ;

    return subset_x ;
  }

  // add two new columns, ids and overlaps
  List out(nc + 2) ;
  CharacterVector onames = df.attr("names") ;
  CharacterVector names(nc + 2) ;

  for (int i = 0; i < nc; i++) {
    out[i] = df[i] ;
    names[i] = onames[i] ;
  }

  // add ids
  out[nc] = ids;
  std::string id_name = ".id_merge" ;
  names[nc] = id_name ;

  // add overlaps
  out[nc + 1] = overlaps ;
  std::string overlaps_name = ".overlap_merge" ;
  names[nc + 1] = overlaps_name;

  out.attr("class") = df.attr("class") ;
  out.attr("row.names") = df.attr("row.names") ;
  out.attr("names") = names ;

  return out ;
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
