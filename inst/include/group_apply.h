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
                         int idx_x, int idx_y) {

  IntegerVector idxs_x = IntegerVector::create(idx_x) ;
  IntegerVector idxs_y = IntegerVector::create(idx_y) ;

  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(idxs_x, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(df_y, names(df_y)).subset(idxs_y, "data.frame");

  int ncols = df_x.size() ;
  bool cols_equal = false;

  for (int i = 0; i < ncols; i++) {
    CharacterVector col_x = subset_x[i] ;
    CharacterVector col_y = subset_y[i] ;

    cols_equal = is_true(all(col_x == col_y)) ;
    if (!cols_equal) break ;
  }
  return cols_equal ;
}


template < typename FN, typename... ARGS >
inline void GroupApply(const GroupedDataFrame& x,
                       const GroupedDataFrame& y,
                       FN&& fn, ARGS&& ... args) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;

  DataFrame labels_x(data_x.attr("labels"));
  DataFrame labels_y(data_y.attr("labels"));

  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for (int nx = 0; nx < ng_x; nx++, ++git_x) {

    SlicingIndex gi_x = *git_x ;

    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for (int ny = 0; ny < ng_y; ny++, ++git_y) {

      SlicingIndex gi_y = *git_y ;

      bool same_groups = compare_rows(labels_x, labels_y, nx, ny);

      if (same_groups) {

        ivl_vector_t vx = makeIntervalVector(data_x, gi_x) ;
        ivl_vector_t vy = makeIntervalVector(data_y, gi_y) ;

        std::bind(std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)...)();
      }
    }
  }
}

#endif
