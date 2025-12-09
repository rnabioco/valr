// dist.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

void dist_grouped(valr::ivl_vector_t& vx, valr::ivl_vector_t& vy, std::vector<int>& indices_x,
                  std::vector<double>& distances, const std::string& dist_fxn) {
  // first build sorted vector of y interval midpoints
  std::vector<int> ref_midpoints;
  for (auto vy_it : vy) {
    int midpoint = (vy_it.start + vy_it.stop) / 2;
    ref_midpoints.push_back(midpoint);
  }

  std::sort(ref_midpoints.begin(), ref_midpoints.end());

  std::size_t low_idx, upper_idx;

  // iterate through x intervals and calculate dist using a binary search
  for (auto vx_it : vx) {
    int midpoint = (vx_it.start + vx_it.stop) / 2;
    auto low_it = std::lower_bound(ref_midpoints.begin(), ref_midpoints.end(), midpoint);

    low_idx = low_it - ref_midpoints.begin();

    if (dist_fxn == "absdist") {
      // set up indexes for closest element, handling edge cases at start and end of x ivl vector
      if (low_idx == 0) {
        // no need to continue return absdist
        int dist = ref_midpoints[low_idx];
        int absdist = std::abs(dist - midpoint);

        distances.push_back(absdist);
        indices_x.push_back(vx_it.value);
        continue;

      } else if (low_idx == ref_midpoints.size()) {
        // just search for closest lower ivl
        low_idx = low_idx - 1;
        upper_idx = low_idx;

      } else {
        // search either
        // get index below and above
        low_idx = low_idx - 1;
        upper_idx = low_idx + 1;
      }

      int left = ref_midpoints[low_idx];
      int right = ref_midpoints[upper_idx];

      int dist_l = std::abs(midpoint - left);
      int dist_r = std::abs(midpoint - right);

      // calc relative distance
      int absdist = std::min(dist_l, dist_r);

      distances.push_back(absdist);
      indices_x.push_back(vx_it.value);

    } else {  // reldist
      // drop intervals at start and end which have no reldist
      if (low_idx == 0 || low_idx == ref_midpoints.size()) {
        continue;
      }

      // get index below and above
      low_idx = low_idx - 1;
      upper_idx = low_idx + 1;

      int left = ref_midpoints[low_idx];
      int right = ref_midpoints[upper_idx];

      int dist_l = std::abs(midpoint - left);
      int dist_r = std::abs(midpoint - right);

      // calc relative distance
      auto reldist = (float)std::min(dist_l, dist_r) / float(right - left);

      distances.push_back(reldist);
      indices_x.push_back(vx_it.value);
    }
  }
}

[[cpp11::register]]
cpp11::writable::data_frame dist_impl(cpp11::data_frame gdf_x, cpp11::data_frame gdf_y,
                                      cpp11::integers x_grp_indexes, cpp11::integers y_grp_indexes,
                                      std::string distcalc) {
  valr::GroupedDataFrame grouped_x(gdf_x);
  valr::GroupedDataFrame grouped_y(gdf_y);

  cpp11::data_frame df_x = grouped_x.data();
  cpp11::data_frame df_y = grouped_y.data();

  std::vector<double> distances;
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

    dist_grouped(vx, vy, indices_x, distances, distcalc);
  }

  // Subset x dataframe by output indices
  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df_x, indices_x);

  // Build output using DataFrameBuilder
  valr::DataFrameBuilder out;

  // Add x dataframe (don't drop chrom)
  out.add_df(subset_x, "", false);

  // Add distances column
  std::string distname = distcalc == "absdist" ? ".absdist" : ".reldist";
  out.add_column(distname, distances);

  int nrows = indices_x.size();
  return out.build(nrows);
}
