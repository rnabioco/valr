// intersect.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void intersect_group(ivl_vector_t vx, ivl_vector_t vy,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& overlap_sizes, bool invert = false) {

  ivl_tree_t tree_y(vy) ;
  ivl_vector_t overlaps ;

  for (auto it : vx) {

    tree_y.findOverlapping(it.start, it.stop, overlaps) ;

    if (overlaps.empty() && invert) {
      indices_x.push_back(it.value) ;
      // store empty placeholer
      indices_y.push_back(NA_INTEGER) ;
      overlap_sizes.push_back(NA_INTEGER) ;
    }

    // store current intervals
    for (auto oit : overlaps) {

      int overlap_size = intervalOverlap(it, oit) ;
      overlap_sizes.push_back(overlap_size) ;

      indices_x.push_back(it.value) ;
      indices_y.push_back(oit.value) ;
    }

    overlaps.clear() ;
  }
}


// [[Rcpp::export]]
DataFrame intersect_impl(GroupedDataFrame x, GroupedDataFrame y,
                         bool invert = false,
                         const std::string& suffix_x = ".x",
                         const std::string& suffix_y = ".y") {

  // indices for subsetting
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;

  // overlap sizes
  std::vector<int> overlap_sizes ;

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  // set up interval trees for each chromosome and apply intersect_group
  GroupApply(x, y, intersect_group, std::ref(indices_x), std::ref(indices_y), std::ref(overlap_sizes), invert);

  DataFrame subset_x = DataFrameSubsetVisitors(data_x, data_x.names()).subset(indices_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(data_y, data_y.names()).subset(indices_y, "data.frame");

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, suffix_x, false) ;

  // y names, data
  out.add_df(subset_y, suffix_y, true) ;

  // overlaps
  out.add_vec(".overlap", wrap(overlap_sizes)) ;

  auto nrows = subset_x.nrows() ;
  auto res = out.format_df(nrows) ;
  return res ;

}

/***R
library(valr)
library(dplyr)

genome <- tibble::tribble(
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
