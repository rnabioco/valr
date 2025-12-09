// gcoverage.cpp
//
// Copyright (C) 2023 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <algorithm>
#include <utility>
#include <vector>

#include "valr/dataframe.hpp"

using namespace cpp11::literals;

typedef std::pair<double, int> posTrack;
typedef std::vector<posTrack> posTracker;

// sort starts and ends while encoding coverage value
// to sum for overlapping intervals
posTracker collatePositions(const std::vector<double>& starts, const std::vector<double>& ends) {
  posTracker p;
  auto n = starts.size();

  if (n != ends.size()) {
    cpp11::stop("incompatible start and end vector supplied");
  }

  for (size_t i = 0; i < n; i++) {
    auto s = starts[i];
    auto e = ends[i];
    p.push_back(posTrack(s, 1));
    p.push_back(posTrack(e, -1));
  }

  std::sort(p.begin(), p.end());

  return p;
}

[[cpp11::register]]
cpp11::writable::data_frame gcoverage_impl(cpp11::data_frame gdf, cpp11::doubles max_coords) {
  valr::GroupedDataFrame grouped(gdf);
  cpp11::data_frame df = grouped.data();

  int ng = grouped.ngroups();
  cpp11::list idx = grouped.indices();

  if (static_cast<int>(max_coords.size()) != ng) {
    cpp11::stop("max_coords must equal the number of groups in data.frame");
  }

  cpp11::doubles df_starts = df["start"];
  cpp11::doubles df_ends = df["end"];

  std::vector<int> out_indices;
  std::vector<int> depths;
  std::vector<double> starts;
  std::vector<double> ends;

  for (int i = 0; i < ng; i++) {
    double max_coord = max_coords[i];

    cpp11::integers indices(idx[i]);
    int ni = indices.size();

    // Extract starts/ends for this group
    std::vector<double> qstarts(ni);
    std::vector<double> qends(ni);

    for (int j = 0; j < ni; j++) {
      int row_idx = indices[j] - 1;  // Convert to 0-based
      qstarts[j] = df_starts[row_idx];
      qends[j] = df_ends[row_idx];
    }

    // track single index of a representative interval
    // this will be used to recover grouped column values
    int grp_idx = indices[0] - 1;

    auto pos = collatePositions(qstarts, qends);

    double prev_pos = 0;
    int cvg = 0;

    for (auto& p : pos) {
      if (p.first > max_coord) {
        cpp11::warning(
            "Out of bounds interval detected at position: %d \n"
            "  Out of bounds intervals will be ignored",
            static_cast<int>(p.first));
        break;
      }

      if (p.first == prev_pos) {
        cvg += p.second;
        continue;
      }
      depths.push_back(cvg);
      starts.push_back(prev_pos);
      ends.push_back(p.first);
      out_indices.push_back(grp_idx);

      prev_pos = p.first;
      cvg += p.second;
    }

    if (ends.back() < max_coord) {
      depths.push_back(0);
      starts.push_back(ends.back());
      ends.push_back(max_coord);
      out_indices.push_back(grp_idx);
    }
  }

  // subset original dataframe, note that the grouped column values will be correct
  // but non-grouped column values are no longer matched to each interval and
  // must be dropped on the R side
  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df, out_indices);

  // Get column names to find start/end positions
  cpp11::strings col_names(subset_x.attr("names"));
  int ncol = subset_x.size();

  // Build output vectors
  int nrow_out = starts.size();
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);
  cpp11::writable::integers new_depths(nrow_out);

  for (int i = 0; i < nrow_out; i++) {
    new_starts[i] = starts[i];
    new_ends[i] = ends[i];
    new_depths[i] = depths[i];
  }

  // Check if .depth already exists in the dataframe
  bool has_depth = false;
  for (int j = 0; j < ncol; j++) {
    if (std::string(col_names[j]) == ".depth") {
      has_depth = true;
      break;
    }
  }

  // Build output list, replacing start/end columns and adding/replacing .depth
  int out_ncol = has_depth ? ncol : ncol + 1;
  cpp11::writable::list out_list(out_ncol);
  cpp11::writable::strings out_names(out_ncol);

  int out_idx = 0;
  for (int j = 0; j < ncol; j++) {
    std::string name(col_names[j]);
    out_names[out_idx] = name;

    if (name == "start") {
      out_list[out_idx] = new_starts;
    } else if (name == "end") {
      out_list[out_idx] = new_ends;
    } else if (name == ".depth") {
      out_list[out_idx] = new_depths;
    } else {
      out_list[out_idx] = subset_x[j];
    }
    out_idx++;
  }

  // Add .depth column if it didn't exist
  if (!has_depth) {
    out_names[out_idx] = ".depth";
    out_list[out_idx] = new_depths;
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}
