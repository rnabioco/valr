// coverage.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void coverage_group(ivl_vector_t vx, ivl_vector_t vy,
                    std::vector<int>& overlap_counts, std::vector<int>& ivls_bases_covered,
                    std::vector<int>& x_ivl_lengths, std::vector<double>& fractions_covered,
                    std::vector<int>& indices_x) {

  ivl_tree_t tree_y(vy) ;
  ivl_vector_t overlaps ;
  IntervalSorterDesc<int, int> intervalSorterDesc;

  for (auto it : vx) {

    indices_x.push_back(it.value);

    tree_y.findOverlapping(it.start, it.stop, overlaps) ;

    // compute number of overlaps
    int overlap_count = overlaps.size();
    overlap_counts.push_back(overlap_count);

    // handle no overlaps and continue
    if (overlap_count == 0) {
      int x_ivl_length = it.stop - it.start ;
      x_ivl_lengths.push_back(x_ivl_length) ;

      ivls_bases_covered.push_back(0) ;
      fractions_covered.push_back(0) ;
      continue;
    }

    // variables to compute number of bases
    int ivl_bases_covered = 0;
    auto x_ivl_start = it.start;
    auto x_ivl_stop = it.stop;

    // total x interval length
    int x_ivl_length = x_ivl_stop - x_ivl_start ;
    x_ivl_lengths.push_back(x_ivl_length) ;
    // merge all overlapping intervals to compute number of bases covered
    // perhaps reuse this code as a seperate function
    // sort overlaps by start descending
    std::sort(overlaps.begin(), overlaps.end(), intervalSorterDesc) ;

    int index = 0; // Stores index of last element

    // Traverse all overlapping intervals
    for (int i = 0; i < overlap_count; i++) {
      // If this is not first Interval and overlaps
      // with the previous one

      if (index != 0 && overlaps[index - 1].start <= overlaps[i].stop)
      {
        while (index != 0 && overlaps[index - 1].start <= overlaps[i].stop)
        {
          // Merge previous and current Intervals
          overlaps[index - 1].stop = std::max(overlaps[index - 1].stop, overlaps[i].stop);
          overlaps[index - 1].start = std::min(overlaps[index - 1].start, overlaps[i].start);
          index--;
        }
      }
      else // Doesn't overlap with previous, add to
        // solution
        overlaps[index] = overlaps[i];

      index++;
    }

    ivl_vector_t mergedOverlaps;
    for (int i = 0; i < index; i++) {
      mergedOverlaps.push_back(overlaps[i]);;
    }
    overlaps.clear();

    // iterate through merged overlaps and compute number of covered bases
    for (auto oit : mergedOverlaps) {

      auto y_ivl_start = oit.start;
      auto y_ivl_stop = oit.stop;

      if (y_ivl_start <  x_ivl_start) {
        y_ivl_start = x_ivl_start ;
      }

      if (y_ivl_stop > x_ivl_stop) {
        y_ivl_stop = x_ivl_stop ;
      }

      int coveredBases = y_ivl_stop - y_ivl_start;
      ivl_bases_covered += coveredBases;

    }

    auto fraction_covered = (double) ivl_bases_covered / x_ivl_length ;
    ivls_bases_covered.push_back(ivl_bases_covered) ;
    fractions_covered.push_back(fraction_covered) ;

    mergedOverlaps.clear();
  }
}


//[[Rcpp::export]]
DataFrame coverage_impl(GroupedDataFrame x, GroupedDataFrame y) {

  // overlapping interval stats
  std::vector<int> overlap_counts ;
  std::vector<int> ivls_bases_covered ;
  std::vector<int> x_ivl_lengths ;
  std::vector<double> fractions_covered ;

  // indices for subsetting
  std::vector<int> indices_x ;

  auto data_x = x.data() ;

  GroupApply(x, y, coverage_group,
             std::ref(overlap_counts), std::ref(ivls_bases_covered),
             std::ref(x_ivl_lengths), std::ref(fractions_covered), std::ref(indices_x));

  // handle condition with empty y df
  // just assign zeros, except for interval length
  if (y.data().nrows() == 0) {
    auto ng_x = x.ngroups() ;

    GroupedDataFrame::group_iterator git_x = x.group_begin() ;
    for (int nx = 0; nx < ng_x; nx++, ++git_x) {

      GroupedSlicingIndex gi_x = *git_x ;
      ivl_vector_t vx = makeIntervalVector(data_x, gi_x) ;

      for (auto it : vx) {
        indices_x.push_back(it.value) ;
        overlap_counts.push_back(0);
        int x_ivl_length = it.stop - it.start ;
        x_ivl_lengths.push_back(x_ivl_length) ;
        ivls_bases_covered.push_back(0) ;
        fractions_covered.push_back(0) ;
      }
    }
  }

  DataFrame subset_x = DataFrameSubsetVisitors(data_x, data_x.names()).subset(indices_x, "data.frame");

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, false) ;

  // additional columns
  out.add_vec(".ints", wrap(overlap_counts)) ;
  out.add_vec(".cov", wrap(ivls_bases_covered)) ;
  out.add_vec(".len", wrap(x_ivl_lengths)) ;
  out.add_vec(".frac", wrap(fractions_covered)) ;

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
