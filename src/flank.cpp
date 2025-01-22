// flank.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void check_coords(double start, double end,
                  double chrom_size, int idx, bool trim,
                  writable::doubles& starts_out,
                  writable::doubles& ends_out,
                  std::vector<int>& df_idx) {

  if (start == end) return ;

  if (start >= 0 && end <= chrom_size) {

    starts_out.push_back(start);
    ends_out.push_back(end);
    df_idx.push_back(idx);

  } else if (trim) {

    if (start < 0) {
      starts_out.push_back(0) ;
    } else {
      starts_out.push_back(start) ;
    }

    if (end > chrom_size) {
      ends_out.push_back(chrom_size) ;
    } else {
      ends_out.push_back(end) ;
    }

    df_idx.push_back(idx);

  } // else trim
}

[[cpp11::register]]
writable::data_frame flank_impl(data_frame df, data_frame genome,
                     double both = 0, double left = 0, double right = 0,
                     bool fraction = false, bool stranded = false, bool trim = false) {

  strings chroms = df["chrom"];
  doubles starts = df["start"];
  doubles ends = df["end"];

  // storage for outputs
  writable::doubles starts_out;
  writable::doubles ends_out;
  std::vector<int> df_idx;

  genome_map_t chrom_sizes = makeChromSizes(genome);
  int lstart, lend, rstart, rend ;

  if (stranded) {

    strings strand = df["strand"];

    for (int i = 0; i < starts.size(); i++) {

      int start = starts[i] ;
      int end = ends[i] ;
      double size = end - start;

      if (fraction) {
        if (strand[i] == "+") {
          lstart = start - std::round(size * left);
          lend = start;
          rstart = end;
          rend = end + std::round(size * right);
        } else {
          lstart = end;
          lend = end + std::round(size * left) ;
          rstart = start - std::round(size * right) ;
          rend = start ;
        }
      } else {
        if (strand[i] == "+") {
          lstart = start - left;
          lend = start;
          rstart = end;
          rend = end + right;
        } else {
          lstart = end;
          lend = end + left;
          rstart = start - right;
          rend = start;
        }
      }

      std::string chrom = chroms[i];
      double chrom_size = chrom_sizes[chrom];

      // check and save coordinates
      check_coords(lstart, lend, chrom_size, i, trim,
                   starts_out, ends_out, df_idx) ;
      check_coords(rstart, rend, chrom_size, i, trim,
                   starts_out, ends_out, df_idx) ;
    }

  } else { // no strand

    for (int i = 0; i < starts.size(); i++) {

      int start = starts[i] ;
      int end = ends[i] ;
      double size = end - start;

      if (fraction) {
        lstart = start - std::round(size * left);
        lend = start;
        rstart = end;
        rend = end + std::round(size * right);
      } else {
        lstart = start - left;
        lend = start;
        rstart = end;
        rend = end + right;
      }

      std::string chrom = chroms[i];
      int chrom_size = chrom_sizes[chrom];

      // check and save coordinates
      check_coords(lstart, lend, chrom_size, i, trim,
                   starts_out, ends_out, df_idx) ;
      check_coords(rstart, rend, chrom_size, i, trim,
                   starts_out, ends_out, df_idx) ;
    }
  }

  writable::data_frame subset = subset_dataframe(df, df_idx) ;

  if (stranded) {
    return writable::data_frame({
      "chrom"_nm = subset["chrom"],
      "start"_nm = starts_out,
      "end"_nm = ends_out,
      "strand"_nm = subset["strand"]
    }) ;
  } else {
    return writable::data_frame({
      "chrom"_nm = subset["chrom"],
      "start"_nm = starts_out,
      "end"_nm = ends_out
    }) ;
  }
}
