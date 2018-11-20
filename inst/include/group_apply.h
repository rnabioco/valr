// group_apply.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__group_apply_H
#define valr__group_apply_H

#include <functional>

#include "valr.h"

inline bool compare_rows(DataFrame df_x, DataFrame df_y,
                         int idx_x, int idx_y, SEXP frame) {

  IntegerVector idxs_x = IntegerVector::create(idx_x) ;
  IntegerVector idxs_y = IntegerVector::create(idx_y) ;

  DataFrame subset_x = subset_dataframe(df_x, idxs_x, frame) ;
  DataFrame subset_y = subset_dataframe(df_y, idxs_y, frame) ;

  // don't compare the .rows column
  int ncols = df_x.size() ;
  CharacterVector cnames_x = df_x.attr("names") ;
  CharacterVector cnames_y = df_y.attr("names") ;
  std::string group_col = ".rows" ;
  bool cols_equal = false;

  for (int i = 0; i < ncols; i++) {
    auto current_col = cnames_x[i];
    if(current_col == group_col) {
      continue ;
    }
    CharacterVector::iterator itr = std::find(cnames_y.begin(), cnames_y.end(), current_col);

    int y_col_idx ;
    if (itr != cnames_y.end()) {
      y_col_idx = std::distance(cnames_y.begin(), itr);
    } else {
      Rcpp::stop("Element not found");
    }

    CharacterVector col_x = subset_x[i] ;
    CharacterVector col_y = subset_y[y_col_idx] ;

    cols_equal = is_true(all(col_x == col_y)) ;
    if (!cols_equal) break ;
  }
  return cols_equal ;
}


template < typename FN, typename... ARGS >
inline void GroupApply(const GroupedDataFrame& x,
                       const GroupedDataFrame& y,
                       SEXP frame,
                       FN&& fn, ARGS&& ... args) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;

  DataFrame labels_x = x.group_data() ;
  DataFrame labels_y = y.group_data() ;

  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for (int nx = 0; nx < ng_x; nx++, ++git_x) {

    GroupedSlicingIndex gi_x = *git_x ;
    if(gi_x.size() == 0) continue ;

    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for (int ny = 0; ny < ng_y; ny++, ++git_y) {

      GroupedSlicingIndex gi_y = *git_y ;
      if(gi_y.size() == 0) continue ;

      bool same_groups = compare_rows(labels_x, labels_y, nx, ny, frame);
      if (same_groups) {

        ivl_vector_t vx = makeIntervalVector(data_x, gi_x) ;
        ivl_vector_t vy = makeIntervalVector(data_y, gi_y) ;

        std::bind(std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)...)();
      }
    }
  }
}

#endif
