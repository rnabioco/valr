// makewindows.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <cmath>
#include <vector>

#include "valr/dataframe.hpp"

using namespace cpp11::literals;

[[cpp11::register]]
cpp11::writable::data_frame makewindows_impl(cpp11::data_frame df, int win_size = 0,
                                             int num_win = 0, int step_size = 0,
                                             bool reverse = false) {
  cpp11::doubles starts = df["start"];
  cpp11::doubles ends = df["end"];

  std::vector<int> starts_out;
  std::vector<int> ends_out;
  std::vector<int> df_idxs;
  std::vector<int> win_ids;

  int nrows = starts.size();

  for (int i = 0; i < nrows; ++i) {
    int start = static_cast<int>(starts[i]);
    int end = static_cast<int>(ends[i]);

    int local_win_size = win_size;
    if (num_win > 0) {
      int ivl_len = end - start;
      if (ivl_len < num_win) {
        continue;
      }
      local_win_size = static_cast<int>(std::round(static_cast<double>(ivl_len) / num_win));
    }

    int by = step_size > 0 ? step_size : local_win_size;

    // Create candidate starts
    std::vector<int> starts_by;
    if (num_win > 0) {
      int win_idx = 0;
      for (int j = start; j < end && win_idx < num_win; j += local_win_size, ++win_idx) {
        starts_by.push_back(j);
      }
    } else {
      for (int j = start; j < end && j + local_win_size - by <= end; j += by) {
        starts_by.push_back(j);
      }
    }

    int nstarts = starts_by.size();
    for (int k = 0; k < nstarts; ++k) {
      int start_by = starts_by[k];
      starts_out.push_back(start_by);

      // For last iteration through starts also extend to end
      if (start_by + local_win_size < end && k < nstarts - 1) {
        ends_out.push_back(start_by + local_win_size);
      } else {
        ends_out.push_back(end);
      }

      if (reverse) {
        win_ids.push_back(nstarts - k);
      } else {
        win_ids.push_back(k + 1);
      }

      df_idxs.push_back(i);
    }
  }

  // Subset input dataframe to get repeated rows
  cpp11::writable::data_frame base_df = valr::subset_dataframe(df, df_idxs);

  // Get column names from base dataframe
  cpp11::strings col_names(base_df.attr("names"));
  int ncol = base_df.size();
  int nrow_out = starts_out.size();

  // Build new vectors for start, end, and .win_id
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);
  cpp11::writable::integers new_win_ids(nrow_out);

  for (int i = 0; i < nrow_out; ++i) {
    new_starts[i] = starts_out[i];
    new_ends[i] = ends_out[i];
    new_win_ids[i] = win_ids[i];
  }

  // Build output list, replacing start/end/.win_id columns
  cpp11::writable::list out_list(ncol);
  cpp11::writable::strings out_names(ncol);

  for (int i = 0; i < ncol; ++i) {
    std::string name(col_names[i]);
    out_names[i] = name;

    if (name == "start") {
      out_list[i] = new_starts;
    } else if (name == "end") {
      out_list[i] = new_ends;
    } else if (name == ".win_id") {
      out_list[i] = new_win_ids;
    } else {
      out_list[i] = base_df[i];
    }
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";

  // Set row names using compact form
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}
