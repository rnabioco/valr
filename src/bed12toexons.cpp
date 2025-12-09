// bed12toexons.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <sstream>
#include <string>
#include <vector>

#include "valr/dataframe.hpp"

using namespace cpp11::literals;

// Parse comma-separated values into a vector of integers
std::vector<int> csv_values(const std::string& csv) {
  std::vector<int> values;
  std::stringstream ss(csv);

  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, ',');

    if (substr.empty())
      break;

    values.push_back(std::atoi(substr.c_str()));
  }

  return values;
}

[[cpp11::register]]
cpp11::writable::data_frame bed12toexons_impl(cpp11::data_frame x) {
  // Input columns (use doubles for coordinates)
  cpp11::doubles starts = x["start"];
  cpp11::strings exon_sizes = x["exon_sizes"];
  cpp11::strings exon_starts = x["exon_starts"];
  cpp11::strings strands = x["strand"];

  // Storage for output
  std::vector<double> starts_out;
  std::vector<double> ends_out;
  std::vector<int> nums_out;  // exon numbers
  std::vector<int> df_idx;

  int nrows = starts.size();

  for (int i = 0; i < nrows; ++i) {
    std::vector<int> exon_start = csv_values(std::string(exon_starts[i]));
    std::vector<int> exon_size = csv_values(std::string(exon_sizes[i]));

    // Calculate starts and ends for each exon
    double start = starts[i];
    int n = exon_start.size();

    for (int j = 0; j < n; ++j) {
      starts_out.push_back(start + exon_start[j]);
      ends_out.push_back(start + exon_start[j] + exon_size[j]);

      if (std::string(strands[i]) == "-") {
        nums_out.push_back(n - j);
      } else {
        nums_out.push_back(j + 1);
      }

      df_idx.push_back(i);
    }
  }

  // Subset input dataframe to get repeated rows
  cpp11::writable::data_frame base_df = valr::subset_dataframe(x, df_idx);

  // Get column names
  cpp11::strings col_names(base_df.attr("names"));
  int ncol = base_df.size();
  int nrow_out = starts_out.size();

  // Build output vectors
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);
  cpp11::writable::integers new_scores(nrow_out);

  for (int i = 0; i < nrow_out; ++i) {
    new_starts[i] = starts_out[i];
    new_ends[i] = ends_out[i];
    new_scores[i] = nums_out[i];
  }

  // Build output list, replacing start/end/score columns
  cpp11::writable::list out_list(ncol);
  cpp11::writable::strings out_names(ncol);

  for (int i = 0; i < ncol; ++i) {
    std::string name(col_names[i]);
    out_names[i] = name;

    if (name == "start") {
      out_list[i] = new_starts;
    } else if (name == "end") {
      out_list[i] = new_ends;
    } else if (name == "score") {
      out_list[i] = new_scores;
    } else {
      out_list[i] = base_df[i];
    }
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}
