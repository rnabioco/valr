// group_apply.h
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__group_apply_H
#define valr__group_apply_H

#include <functional>

#include "valr.h"

template < typename FN, typename... ARGS >
inline void GroupApply(const ValrGroupedDataFrame& x,
                       const ValrGroupedDataFrame& y,
                       const IntegerVector& shared_grps_x,
                       const IntegerVector& shared_grps_y,
                       FN&& fn, ARGS&& ... args) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  int ng_x = shared_grps_x.size() ;
  int ng_y = shared_grps_y.size() ;

  if (ng_x != ng_y) {
    stop("incompatible groups found between x and y dataframes") ;
  }

  // access the group .rows list
  ListView grp_indices_x(x.indices()) ;
  ListView grp_indices_y(y.indices()) ;

  for (int i = 0; i < ng_x; i++) {
    // get next row index to subset from x and y groups
    // convert from R to Cpp index
    int shared_x_index = shared_grps_x[i] - 1;
    int shared_y_index = shared_grps_y[i] - 1;

    // subset the group lists and return integervectors
    IntegerVector gi_x, gi_y ;
    gi_x = grp_indices_x[shared_x_index] ;
    gi_y = grp_indices_y[shared_y_index] ;

    if (gi_x.size() == 0 || gi_y.size() == 0) {
      continue ;
    }

    ivl_vector_t vx = makeIntervalVector(data_x, gi_x) ;
    ivl_vector_t vy = makeIntervalVector(data_y, gi_y) ;

    std::bind(std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)...)();
  }
}


#endif
