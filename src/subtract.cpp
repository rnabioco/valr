// subtract.cpp
//
// Copyright (C) 2016 - 2022 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void subtract_group(ivl_vector_t vx, ivl_vector_t vy,
                    std::vector<int>& indices_out,
                    std::vector<int>& starts_out, std::vector<int>& ends_out) {

  ivl_tree_t tree_y(std::move(vy)) ;
  ivl_vector_t overlaps ;
  IntervalStartSorter<int, int> ivl_sorter ;

  for (auto it : vx) {

    auto x_start = it.start;
    auto x_stop = it.stop;

    overlaps = tree_y.findOverlapping(it.start, it.stop) ;

    // compute number of overlaps
    int overlap_count = overlaps.size();

    // handle no overlaps and continue
    if (overlap_count == 0) {
      indices_out.push_back(it.value) ;
      starts_out.push_back(it.start) ;
      ends_out.push_back(it.stop) ;
      continue;
    }

    // sort overlaps, as sort order not guaranteed
    std::sort(overlaps.begin(), overlaps.end(), ivl_sorter) ;

    // iterate through overlaps with current x  interval
    // modifying start and stop as necessary
    for (auto oit : overlaps) {

      auto y_start = oit.start;
      auto y_stop = oit.stop;

      if (x_start > x_stop) {
        break ;
      } else if (y_start <= x_start) {
        // advance x_start to end of y unless y is shorter than x
        if (x_start > y_stop) {
          continue ;
        } else {
          x_start = y_stop ;
        };
      } else {
        // report new interval
        indices_out.push_back(it.value) ;
        starts_out.push_back(x_start) ;
        ends_out.push_back(y_start) ;
        // advance to end of y ivl
        x_start = y_stop;
      }
    }

    if (x_start < x_stop) {
      // report interval
      indices_out.push_back(it.value) ;
      starts_out.push_back(x_start) ;
      ends_out.push_back(x_stop) ;
    }

    overlaps.clear() ;
  }
}

//[[Rcpp::export]]
DataFrame subtract_impl(ValrGroupedDataFrame gdf_x, ValrGroupedDataFrame gdf_y,
                        IntegerVector x_grp_indexes,
                        IntegerVector y_grp_indexes) {

  std::vector<std::string> chrom_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;

  DataFrame df_x = gdf_x.data() ;
  DataFrame df_y = gdf_y.data() ;

  // indices_to_report
  std::vector<int> indices_out ;

  // set up interval trees for each chromosome and apply subtract_group
  GroupApply(gdf_x, gdf_y, x_grp_indexes, y_grp_indexes, subtract_group, std::ref(indices_out), std::ref(starts_out), std::ref(ends_out));

  // extract out x data, new intervals will be generated as copies of the parent interval
  DataFrame out = subset_dataframe(df_x, indices_out) ;

  // assign new starts and ends
  out["start"] = starts_out ;
  out["end"] = ends_out ;

  return out ;

}

/***R
library(dplyr)
library(valr)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7
)

n <- 1e4
x <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)
y <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)

subtract_impl(x, y) %>% as_data_frame()

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    200
) %>% group_by(chrom)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1000,    2000
) %>% group_by(chrom)

subtract_impl(x, y)

*/
