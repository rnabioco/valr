// coverage.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <algorithm>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

void coverage_group(valr::ivl_vector_t vx, valr::ivl_vector_t vy, std::vector<int>& overlap_counts,
                    std::vector<int>& ivls_bases_covered, std::vector<int>& x_ivl_lengths,
                    std::vector<double>& fractions_covered, std::vector<int>& indices_x,
                    int min_overlap) {
  valr::ivl_tree_t tree_y(std::move(vy));
  valr::ivl_vector_t overlaps;
  valr::IntervalStopDescCmp<int, int> intervalSorterDesc;

  for (const auto& it : vx) {
    indices_x.push_back(it.value);

    overlaps = tree_y.findOverlapping(it.start, it.stop, min_overlap);

    // compute number of overlaps
    int overlap_count = overlaps.size();
    overlap_counts.push_back(overlap_count);

    // handle no overlaps and continue
    if (overlap_count == 0) {
      int x_ivl_length = it.stop - it.start;
      x_ivl_lengths.push_back(x_ivl_length);

      ivls_bases_covered.push_back(0);
      fractions_covered.push_back(0);
      continue;
    }

    // variables to compute number of bases
    int ivl_bases_covered = 0;
    auto x_ivl_start = it.start;
    auto x_ivl_stop = it.stop;

    // total x interval length
    int x_ivl_length = x_ivl_stop - x_ivl_start;
    x_ivl_lengths.push_back(x_ivl_length);

    // merge all overlapping intervals to compute number of bases covered
    // sort overlaps by start descending
    std::sort(overlaps.begin(), overlaps.end(), intervalSorterDesc);

    int index = 0;  // Stores index of last element

    // Traverse all overlapping intervals
    for (int i = 0; i < overlap_count; i++) {
      // If this is not first Interval and overlaps
      // with the previous one

      if (index != 0 && overlaps[index - 1].start <= overlaps[i].stop) {
        while (index != 0 && overlaps[index - 1].start <= overlaps[i].stop) {
          // Merge previous and current Intervals
          overlaps[index - 1].stop = std::max(overlaps[index - 1].stop, overlaps[i].stop);
          overlaps[index - 1].start = std::min(overlaps[index - 1].start, overlaps[i].start);
          index--;
        }
      } else {
        // Doesn't overlap with previous, add to solution
        overlaps[index] = overlaps[i];
      }

      index++;
    }

    valr::ivl_vector_t mergedOverlaps;
    for (int i = 0; i < index; i++) {
      mergedOverlaps.push_back(overlaps[i]);
    }
    overlaps.clear();

    // iterate through merged overlaps and compute number of covered bases
    for (const auto& oit : mergedOverlaps) {
      auto y_ivl_start = oit.start;
      auto y_ivl_stop = oit.stop;

      if (y_ivl_start < x_ivl_start) {
        y_ivl_start = x_ivl_start;
      }

      if (y_ivl_stop > x_ivl_stop) {
        y_ivl_stop = x_ivl_stop;
      }

      int coveredBases = y_ivl_stop - y_ivl_start;
      ivl_bases_covered += coveredBases;
    }

    auto fraction_covered = (double)ivl_bases_covered / x_ivl_length;
    ivls_bases_covered.push_back(ivl_bases_covered);
    fractions_covered.push_back(fraction_covered);

    mergedOverlaps.clear();
  }
}

[[cpp11::register]]
cpp11::writable::data_frame coverage_impl(cpp11::data_frame gdf_x, cpp11::data_frame gdf_y,
                                          cpp11::integers x_grp_indexes,
                                          cpp11::integers y_grp_indexes, int min_overlap) {
  valr::GroupedDataFrame grouped_x(gdf_x);
  valr::GroupedDataFrame grouped_y(gdf_y);

  cpp11::data_frame df_x = grouped_x.data();
  cpp11::data_frame df_y = grouped_y.data();

  // overlapping interval stats
  std::vector<int> overlap_counts;
  std::vector<int> ivls_bases_covered;
  std::vector<int> x_ivl_lengths;
  std::vector<double> fractions_covered;

  // indices for subsetting
  std::vector<int> indices_x;

  // Get group indices lists
  cpp11::list idx_x = grouped_x.indices();
  cpp11::list idx_y = grouped_y.indices();

  int ng = x_grp_indexes.size();

  // Iterate over shared groups
  for (int i = 0; i < ng; i++) {
    // Convert from R to C++ indexing
    int shared_x_index = x_grp_indexes[i] - 1;
    int shared_y_index = y_grp_indexes[i] - 1;

    cpp11::integers gi_x(idx_x[shared_x_index]);
    cpp11::integers gi_y(idx_y[shared_y_index]);

    if (gi_x.size() == 0 || gi_y.size() == 0) {
      continue;
    }

    valr::ivl_vector_t vx = valr::makeIntervalVector(df_x, gi_x);
    valr::ivl_vector_t vy = valr::makeIntervalVector(df_y, gi_y);

    coverage_group(vx, vy, overlap_counts, ivls_bases_covered, x_ivl_lengths, fractions_covered,
                   indices_x, min_overlap);
  }

  // handle condition with empty y df
  // just assign zeros, except for interval length
  if (df_y.nrow() == 0) {
    int ng_x = grouped_x.ngroups();

    for (int nx = 0; nx < ng_x; nx++) {
      cpp11::integers gi_x(idx_x[nx]);
      valr::ivl_vector_t vx = valr::makeIntervalVector(df_x, gi_x);

      for (const auto& it : vx) {
        indices_x.push_back(it.value);
        overlap_counts.push_back(0);
        int x_ivl_length = it.stop - it.start;
        x_ivl_lengths.push_back(x_ivl_length);
        ivls_bases_covered.push_back(0);
        fractions_covered.push_back(0);
      }
    }
  }

  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df_x, indices_x);

  valr::DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, "", false);

  // additional columns
  out.names.push_back(".ints");
  cpp11::writable::integers ints_vec(overlap_counts.size());
  for (size_t i = 0; i < overlap_counts.size(); ++i) {
    ints_vec[i] = overlap_counts[i];
  }
  out.data.push_back(ints_vec);

  out.names.push_back(".cov");
  cpp11::writable::integers cov_vec(ivls_bases_covered.size());
  for (size_t i = 0; i < ivls_bases_covered.size(); ++i) {
    cov_vec[i] = ivls_bases_covered[i];
  }
  out.data.push_back(cov_vec);

  out.names.push_back(".len");
  cpp11::writable::integers len_vec(x_ivl_lengths.size());
  for (size_t i = 0; i < x_ivl_lengths.size(); ++i) {
    len_vec[i] = x_ivl_lengths[i];
  }
  out.data.push_back(len_vec);

  out.names.push_back(".frac");
  cpp11::writable::doubles frac_vec(fractions_covered.size());
  for (size_t i = 0; i < fractions_covered.size(); ++i) {
    frac_vec[i] = fractions_covered[i];
  }
  out.data.push_back(frac_vec);

  int nrows = indices_x.size();
  return out.build(nrows);
}
