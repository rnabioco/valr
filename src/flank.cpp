// flank.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/genome.hpp"

using namespace cpp11::literals;

void check_coords(int start, int end, int chrom_size, int idx, bool trim,
                  std::vector<int>& starts_out, std::vector<int>& ends_out,
                  std::vector<int>& df_idx) {
  if (start == end)
    return;

  if (start >= 0 && end <= chrom_size) {
    starts_out.push_back(start);
    ends_out.push_back(end);
    df_idx.push_back(idx);

  } else if (trim) {
    if (start < 0) {
      starts_out.push_back(0);
    } else {
      starts_out.push_back(start);
    }

    if (end > chrom_size) {
      ends_out.push_back(chrom_size);
    } else {
      ends_out.push_back(end);
    }

    df_idx.push_back(idx);
  }
}

[[cpp11::register]]
cpp11::writable::data_frame flank_impl(cpp11::data_frame df, cpp11::data_frame genome,
                                       double both = 0, double left = 0, double right = 0,
                                       bool fraction = false, bool stranded = false,
                                       bool trim = false) {
  cpp11::strings chroms = df["chrom"];
  cpp11::doubles starts = df["start"];
  cpp11::doubles ends = df["end"];

  // Storage for outputs
  std::vector<int> starts_out;
  std::vector<int> ends_out;
  std::vector<int> df_idx;

  valr::genome_map_t chrom_sizes = valr::make_genome_map(genome);
  int lstart, lend, rstart, rend;

  int nrows = starts.size();

  if (stranded) {
    cpp11::strings strands = df["strand"];

    for (int i = 0; i < nrows; ++i) {
      int start = starts[i];
      int end = ends[i];
      double size = end - start;

      if (fraction) {
        if (std::string(strands[i]) == "+") {
          lstart = start - static_cast<int>(std::round(size * left));
          lend = start;
          rstart = end;
          rend = end + static_cast<int>(std::round(size * right));
        } else {
          lstart = end;
          lend = end + static_cast<int>(std::round(size * left));
          rstart = start - static_cast<int>(std::round(size * right));
          rend = start;
        }
      } else {
        if (std::string(strands[i]) == "+") {
          lstart = start - static_cast<int>(left);
          lend = start;
          rstart = end;
          rend = end + static_cast<int>(right);
        } else {
          lstart = end;
          lend = end + static_cast<int>(left);
          rstart = start - static_cast<int>(right);
          rend = start;
        }
      }

      std::string chrom(chroms[i]);
      int chrom_size = chrom_sizes[chrom];

      // Check and save coordinates
      check_coords(lstart, lend, chrom_size, i, trim, starts_out, ends_out, df_idx);
      check_coords(rstart, rend, chrom_size, i, trim, starts_out, ends_out, df_idx);
    }

  } else {  // no strand
    for (int i = 0; i < nrows; ++i) {
      int start = starts[i];
      int end = ends[i];
      double size = end - start;

      if (fraction) {
        lstart = start - static_cast<int>(std::round(size * left));
        lend = start;
        rstart = end;
        rend = end + static_cast<int>(std::round(size * right));
      } else {
        lstart = start - static_cast<int>(left);
        lend = start;
        rstart = end;
        rend = end + static_cast<int>(right);
      }

      std::string chrom(chroms[i]);
      int chrom_size = chrom_sizes[chrom];

      // Check and save coordinates
      check_coords(lstart, lend, chrom_size, i, trim, starts_out, ends_out, df_idx);
      check_coords(rstart, rend, chrom_size, i, trim, starts_out, ends_out, df_idx);
    }
  }

  // Subset input dataframe
  cpp11::writable::data_frame base_df = valr::subset_dataframe(df, df_idx);

  // Get column names
  cpp11::strings col_names(base_df.attr("names"));
  int ncol = base_df.size();
  int nrow_out = starts_out.size();

  // Build output vectors (as doubles to match input type)
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);

  for (int i = 0; i < nrow_out; ++i) {
    new_starts[i] = starts_out[i];
    new_ends[i] = ends_out[i];
  }

  // Build output list, replacing start/end columns
  cpp11::writable::list out_list(ncol);
  cpp11::writable::strings out_names(ncol);

  for (int i = 0; i < ncol; ++i) {
    std::string name(col_names[i]);
    out_names[i] = name;

    if (name == "start") {
      out_list[i] = new_starts;
    } else if (name == "end") {
      out_list[i] = new_ends;
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
