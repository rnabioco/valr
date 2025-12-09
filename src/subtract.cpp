// subtract.cpp
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

// Process subtraction for a single group
void subtract_group(valr::ivl_vector_t vx, valr::ivl_vector_t vy, std::vector<int>& indices_out,
                    std::vector<double>& starts_out, std::vector<double>& ends_out,
                    int min_overlap) {
  valr::ivl_tree_t tree_y(std::move(vy));
  valr::ivl_vector_t overlaps;

  for (const auto& it : vx) {
    auto x_start = it.start;
    auto x_stop = it.stop;

    overlaps = tree_y.findOverlapping(it.start, it.stop, min_overlap);

    // compute number of overlaps
    int overlap_count = overlaps.size();

    // handle no overlaps and continue
    if (overlap_count == 0) {
      indices_out.push_back(it.value);
      starts_out.push_back(it.start);
      ends_out.push_back(it.stop);
      continue;
    }

    // sort overlaps by start position
    std::sort(overlaps.begin(), overlaps.end(), valr::IntervalStartCmp<int, int>());

    // iterate through overlaps with current x interval
    // modifying start and stop as necessary
    for (const auto& oit : overlaps) {
      auto y_start = oit.start;
      auto y_stop = oit.stop;

      if (x_start > x_stop) {
        break;
      } else if (y_start <= x_start) {
        // advance x_start to end of y unless y is shorter than x
        if (x_start > y_stop) {
          continue;
        } else {
          x_start = y_stop;
        }
      } else {
        // report new interval
        indices_out.push_back(it.value);
        starts_out.push_back(x_start);
        ends_out.push_back(y_start);
        // advance to end of y ivl
        x_start = y_stop;
      }
    }

    if (x_start < x_stop) {
      // report interval
      indices_out.push_back(it.value);
      starts_out.push_back(x_start);
      ends_out.push_back(x_stop);
    }

    overlaps.clear();
  }
}

[[cpp11::register]]
cpp11::writable::data_frame subtract_impl(cpp11::data_frame gdf_x, cpp11::data_frame gdf_y,
                                          cpp11::integers x_grp_indexes,
                                          cpp11::integers y_grp_indexes, int min_overlap) {
  valr::GroupedDataFrame grouped_x(gdf_x);
  valr::GroupedDataFrame grouped_y(gdf_y);

  cpp11::data_frame df_x = grouped_x.data();
  cpp11::data_frame df_y = grouped_y.data();

  // Output vectors
  std::vector<int> indices_out;
  std::vector<double> starts_out;
  std::vector<double> ends_out;

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

    subtract_group(vx, vy, indices_out, starts_out, ends_out, min_overlap);
  }

  // Subset x dataframe by output indices
  cpp11::writable::data_frame out = valr::subset_dataframe(df_x, indices_out);

  // Get column names to find start/end positions
  cpp11::strings col_names(out.attr("names"));
  int ncol = out.size();
  int nrow_out = starts_out.size();

  // Build new start/end vectors
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);

  for (int i = 0; i < nrow_out; i++) {
    new_starts[i] = starts_out[i];
    new_ends[i] = ends_out[i];
  }

  // Build output list, replacing start/end columns
  cpp11::writable::list out_list(ncol);
  cpp11::writable::strings out_names(ncol);

  for (int j = 0; j < ncol; j++) {
    std::string name(col_names[j]);
    out_names[j] = name;

    if (name == "start") {
      out_list[j] = new_starts;
    } else if (name == "end") {
      out_list[j] = new_ends;
    } else {
      out_list[j] = out[j];
    }
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}
