// intervals.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__intervals_H
#define valr__intervals_H

#include "valr.h"

// main interval types used in valr
typedef Interval<int>      ivl_t ;
typedef std::vector<ivl_t> ivl_vector_t ;
typedef IntervalTree<int>  ivl_tree_t ;

// the value field of intervals in the returned vector correspond to the index
// of the interval in the original dataframe (i.e., the values of the
// SlicingIndex)

inline ivl_vector_t makeIntervalVector(DataFrame df, GroupedSlicingIndex si,
                                       std::string col_start = "start",
                                       std::string col_end = "end") {

  ivl_vector_t ivls ;

  IntegerVector starts = df[col_start] ;
  IntegerVector ends   = df[col_end] ;

  int size = si.size() ;

  for (int i = 0; i < size; ++i) {
    int j = si[i] ;
    ivls.push_back(ivl_t(starts[j], ends[j], j)) ;
  }
  return ivls ;
}

#endif
