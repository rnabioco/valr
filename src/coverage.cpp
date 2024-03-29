// coverage.cpp
//
// Copyright (C) 2016 - 2022 Jay Hesselberth and Kent Riemondy
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

  ivl_tree_t tree_y(std::move(vy)) ;
  ivl_vector_t overlaps ;
  IntervalSorterDesc<int, int> intervalSorterDesc;

  for (auto it : vx) {

    indices_x.push_back(it.value);

    overlaps = tree_y.findOverlapping(it.start, it.stop) ;

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
DataFrame coverage_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y,
                        IntegerVector x_grp_indexes,
                        IntegerVector y_grp_indexes) {

  // overlapping interval stats
  std::vector<int> overlap_counts ;
  std::vector<int> ivls_bases_covered ;
  std::vector<int> x_ivl_lengths ;
  std::vector<double> fractions_covered ;

  // indices for subsetting
  std::vector<int> indices_x ;

  auto data_x = x.data() ;

  GroupApply(x, y, x_grp_indexes, y_grp_indexes,
             coverage_group,
             std::ref(overlap_counts), std::ref(ivls_bases_covered),
             std::ref(x_ivl_lengths), std::ref(fractions_covered), std::ref(indices_x));

  // handle condition with empty y df
  // just assign zeros, except for interval length
  if (y.data().nrows() == 0) {
    auto ng_x = x.ngroups() ;

    ListView idx_x(x.indices()) ;

    for (int nx = 0; nx < ng_x; nx++) {
      IntegerVector gi_x ;
      gi_x = idx_x[nx];
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

  DataFrame subset_x = subset_dataframe(data_x, indices_x) ;

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, false) ;

  // additional columns
  out.names.push_back(".ints") ;
  out.data.push_back(overlap_counts) ;

  out.names.push_back(".cov") ;
  out.data.push_back(ivls_bases_covered) ;

  out.names.push_back(".len") ;
  out.data.push_back(x_ivl_lengths) ;

  out.names.push_back(".frac") ;
  out.data.push_back(fractions_covered) ;

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
