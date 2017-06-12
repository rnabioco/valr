// flank.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

void check_coords(int start, int end,
                  int chrom_size, int idx, bool trim,
                  std::vector<int>& starts_out,
                  std::vector<int>& ends_out,
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

//[[Rcpp::export]]
DataFrame flank_impl(DataFrame df, DataFrame genome,
                     double both = 0, double left = 0, double right = 0,
                     bool fraction = false, bool stranded = false, bool trim = false) {

  std::vector<std::string> chroms = df["chrom"];
  IntegerVector starts = df["start"];
  IntegerVector ends = df["end"];

  // storage for outputs
  std::vector<int> starts_out;
  std::vector<int> ends_out;
  std::vector<int> df_idx;

  genome_map_t chrom_sizes = makeChromSizes(genome);
  int lstart, lend, rstart, rend ;

  if (stranded) {

    std::vector<std::string> strands = df["strand"];

    for (int i = 0; i < starts.size(); i++) {

      int start = starts[i] ;
      int end = ends[i] ;
      double size = end - start;

      if (fraction) {
        if (strands[i] == "+") {
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
        if (strands[i] == "+") {
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
      int chrom_size = chrom_sizes[chrom];

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

  DataFrame out = DataFrameSubsetVisitors(df, names(df)).subset(df_idx, "data.frame");

  out["start"] = starts_out;
  out["end"] = ends_out;

  return out;
}


/*** R
library(valr)
library(dplyr)

genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))
x <- bed_random(genome)

devtools::load_all()
flank_impl(x, genome, both = 100) %>% as_data_frame()
*/

