// makewindows.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

std::vector<int> seq_by(int from, int to, int by, int win_size) {

  // determine output size
  std::size_t n = ((to - from) / by) + 1 ;

  // generate output vector
  IntegerVector res = Rcpp::rep(from, n) ;
  std::vector<int> out ;

  // iterate through from to to by step size
  int idx = 0 ;
  for (int i = from; i < to && i + win_size <= to; i += by) {
    out.push_back(i) ;
  }

  return out ;
}

//[[Rcpp::export]]
DataFrame makewindows_impl(DataFrame df, int win_size = 0, int num_win = 0,
                           int step_size = 0, bool reverse = false) {

  NumericVector starts = df["start"] ;
  NumericVector ends = df["end"] ;

  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
  std::vector<int> df_idxs ;
  std::vector<int> win_ids;

  // used in loops
  std::vector<int> ends_by ;
  std::vector<int> ids ;

  for (int i = 0; i < starts.size(); ++i) {

    auto start = starts[i] ;
    auto end = ends[i] ;

    if (num_win > 0) {
      win_size = round((end - start) / num_win) ;
    }

    int by = win_size - step_size ;

    // create candidate starts
    std::vector<int> starts_by = seq_by(start, end, by, win_size) ;

    for (int j = 0; j < starts_by.size(); ++j) {

      auto start_by = starts_by[j] ;
      auto end_by = start + win_size ;

      if (end_by < end) {

        if (end_by + win_size > end) {
          ends_by.push_back(end) ;
        } else {
          ends_by.push_back(start_by + win_size) ;
        }

        ids.push_back(j + 1) ;
      }
    }

    if (reverse) {
      std::reverse(ids.begin(), ids.end());
    }

    // nums of reps of current x idx
    IntegerVector x_idxs = Rcpp::rep(i, starts_by.size()) ;

    // add new ivls
    starts_out.insert(starts_out.end(), starts_by.begin(), starts_by.end());
    ends_out.insert(ends_out.end(), ends_by.begin(), ends_by.end());
    win_ids.insert(win_ids.end(), ids.begin(), ids.end());

    df_idxs.insert(df_idxs.end(), x_idxs.begin(), x_idxs.end());

    // reset
    ends_by.clear() ;
    ids.clear() ;
  }

  DataFrame out = DataFrameSubsetVisitors(df, names(df)).subset(df_idxs, "data.frame");

  // replace original starts, ends, and .win_id
  out["start"] = starts_out ;
  out["end"] = ends_out ;
  out[".win_id"] = win_ids ;

  return out ;
}

/*** R
library(dplyr)
x <- trbl_interval(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 100,    200,  'A',   '.',    '+'
)

devtools::load_all()
makewindows_impl(x, win_size = 10) %>% as_data_frame()
*/
