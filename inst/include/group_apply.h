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

std::vector<int> shared_row_indexes(const ValrGroupedDataFrame& x,
                                    const ValrGroupedDataFrame& y) ;

DataFrame rowwise_subset_df(const DataFrame& x,
                            IntegerVector row_indices) ;

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
inline void GroupApply(const ValrGroupedDataFrame& x,
                       const ValrGroupedDataFrame& y,
                       IntegerVector grp_idx_x,
                       IntegerVector grp_idx_y,
                       FN&& fn, ARGS&& ... args) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  if(grp_idx_x.size() != grp_idx_y.size()){
    stop("incompatible groups found between x and y dataframes") ;
  }

  ListView idx_x(x.indices()) ;
  ListView idx_y(y.indices()) ;

  int ng_x = grp_idx_x.size() ;
  for (int i = 0; i < ng_x; i++) {
    int x_idx = grp_idx_x[i] ;
    int y_idx = grp_idx_y[i] ;

    IntegerVector gi_x, gi_y ;
    gi_x = idx_x[x_idx - 1];
    gi_y = idx_y[y_idx - 1] ;

    if(gi_x.size() == 0 || gi_y.size() == 0) {
      continue ;
    }

    ivl_vector_t vx = makeIntervalVector(data_x, gi_x) ;
    ivl_vector_t vy = makeIntervalVector(data_y, gi_y) ;

    std::bind(std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)...)();
  }
}


#endif
