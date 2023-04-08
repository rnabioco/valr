// grouped_dataframe.h
//
// Copyright (C) 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__grouped_dataframe_H
#define valr__grouped_dataframe_H

#include "valr.h"

class ValrGroupedDataFrame;

class ValrGroupedDataFrame {
public:

  ValrGroupedDataFrame(DataFrame x) ;

  DataFrame& data() {
    return data_;
  }
  const DataFrame& data() const {
    return data_;
  }

  inline int size() const {
    return data_.size();
  }

  inline int ngroups() const {
    return groups.nrows();
  }

  inline int nrows() const {
    return data_.nrows();
  }

  inline const DataFrame& group_data() const {
    return groups ;
  }

  inline SEXP indices() const {
    return groups[groups.size() - 1] ;
  }

  template <typename Data>
  static void strip_groups(Data& x) {
    x.attr("groups") = R_NilValue;
  }

private:
  DataFrame data_;
  DataFrame groups;
};

DataFrame extract_groups(const DataFrame& x) ;

DataFrame rowwise_subset_df(const DataFrame& x,
                            IntegerVector row_indices) ;

DataFrame rowwise_subset_df(const DataFrame& x,
                            std::vector<int> row_indices) ;

inline bool compare_rows(DataFrame df_x, DataFrame df_y,
                         int idx_x, int idx_y) {

  IntegerVector idxs_x = IntegerVector::create(idx_x) ;
  IntegerVector idxs_y = IntegerVector::create(idx_y) ;

  DataFrame subset_x = subset_dataframe(df_x, idxs_x) ;
  DataFrame subset_y = subset_dataframe(df_y, idxs_y) ;

  // don't compare the .rows column
  int ncols = df_x.size() ;
  CharacterVector cnames_x = df_x.attr("names") ;
  CharacterVector cnames_y = df_y.attr("names") ;
  std::string group_col = ".rows" ;
  bool cols_equal = false;

  for (int i = 0; i < ncols; i++) {
    auto current_col = cnames_x[i];
    if (current_col == group_col) {
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

#endif

