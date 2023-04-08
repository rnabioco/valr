// IntervalTree_ext.h
//
// Copyright (C) 2022 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.
//
// intervalOverlap and IntervalStartSorter
// from https://github.com/ekg/intervaltree
// circa commit 8fc4be9 Dec 13, 2015
// see inst/include/intervalTree.h for copyright
//

#ifndef valr__intervalTree_ext_H
#define valr__intervalTree_ext_H

#include "valr.h"

template <class T, typename K = int>
K intervalOverlap(const Interval<T, K>& a, const Interval<T, K>& b) {
  return std::min(a.stop, b.stop) - std::max(a.start, b.start) ;
}

template <class T, typename K = int>
class IntervalStartSorter {
public:
  bool operator()(const Interval<T, K>& a, const Interval<T, K>& b) {
    return a.start < b.start;
  }
};

template <class T, typename K = int>
class IntervalSorterDesc {
public:
  bool operator()(const Interval<T, K>& a, const Interval<T, K>& b) {
    return a.start > b.start;
  }
};

#endif
