// intersect.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <string>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

void intersect_group(valr::ivl_vector_t vx, valr::ivl_vector_t vy, std::vector<int>& indices_x,
                     std::vector<int>& indices_y, std::vector<int>& overlap_sizes, bool invert,
                     int min_overlap) {
  // do not call vy after std::move
  valr::ivl_tree_t tree_y(std::move(vy));

  for (const auto& it : vx) {
    bool found_overlap = false;

    // Use visitor pattern to avoid allocating a new vector per query
    tree_y.visit_overlapping(
        it.start, it.stop,
        [&](const valr::ivl_t& oit) {
          found_overlap = true;
          int overlap_size = valr::intervalOverlap(it, oit);
          overlap_sizes.push_back(overlap_size);
          indices_x.push_back(it.value);
          indices_y.push_back(oit.value);
        },
        min_overlap);

    if (!found_overlap && invert) {
      indices_x.push_back(it.value);
      indices_y.push_back(NA_INTEGER);
      overlap_sizes.push_back(NA_INTEGER);
    }
  }
}

// Find x intervals in groups not present in y
void unmatched_groups(valr::GroupedDataFrame& grouped_x, valr::GroupedDataFrame& grouped_y,
                      std::vector<int>& indices_x, std::vector<int>& indices_y,
                      std::vector<int>& overlap_sizes) {
  int ng_x = grouped_x.ngroups();
  int ng_y = grouped_y.ngroups();

  cpp11::data_frame labels_x = grouped_x.group_data();
  cpp11::data_frame labels_y = grouped_y.group_data();

  cpp11::list idx_x = grouped_x.indices();

  for (int nx = 0; nx < ng_x; nx++) {
    cpp11::integers gi_x(idx_x[nx]);

    bool match = false;

    for (int ny = 0; ny < ng_y; ny++) {
      match = valr::compare_rows(labels_x, labels_y, nx, ny);
      if (match)
        break;
    }

    if (match)
      continue;

    for (int i = 0; i < gi_x.size(); i++) {
      indices_x.push_back(gi_x[i] - 1);
      indices_y.push_back(NA_INTEGER);
      overlap_sizes.push_back(NA_INTEGER);
    }
  }
}

[[cpp11::register]]
cpp11::writable::data_frame intersect_impl(cpp11::data_frame gdf_x, cpp11::data_frame gdf_y,
                                           cpp11::integers x_grp_indexes,
                                           cpp11::integers y_grp_indexes, bool invert,
                                           std::string suffix_x, std::string suffix_y,
                                           int min_overlap) {
  valr::GroupedDataFrame grouped_x(gdf_x);
  valr::GroupedDataFrame grouped_y(gdf_y);

  cpp11::data_frame df_x = grouped_x.data();
  cpp11::data_frame df_y = grouped_y.data();

  // indices for subsetting
  std::vector<int> indices_x;
  std::vector<int> indices_y;

  // overlap sizes
  std::vector<int> overlap_sizes;

  // Pre-reserve capacity to reduce reallocations
  // Estimate ~2 overlaps per x interval on average
  size_t estimated_size = static_cast<size_t>(df_x.nrow()) * 2;
  indices_x.reserve(estimated_size);
  indices_y.reserve(estimated_size);
  overlap_sizes.reserve(estimated_size);

  // find unmatched x intervals
  if (invert) {
    unmatched_groups(grouped_x, grouped_y, indices_x, indices_y, overlap_sizes);
  }

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

    intersect_group(vx, vy, indices_x, indices_y, overlap_sizes, invert, min_overlap);
  }

  // Subset dataframes by output indices
  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df_x, indices_x);
  cpp11::writable::data_frame subset_y = valr::subset_dataframe(df_y, indices_y);

  // Build output using DataFrameBuilder
  valr::DataFrameBuilder out;

  // x names, data (don't drop chrom)
  out.add_df(subset_x, suffix_x, false);

  // y names, data (drop chrom)
  out.add_df(subset_y, suffix_y, true);

  // overlaps
  out.names.push_back(".overlap");
  cpp11::writable::integers overlap_vec(overlap_sizes.size());
  for (size_t i = 0; i < overlap_sizes.size(); ++i) {
    overlap_vec[i] = overlap_sizes[i];
  }
  out.data.push_back(overlap_vec);

  int nrows = indices_x.size();
  return out.build(nrows);
}
