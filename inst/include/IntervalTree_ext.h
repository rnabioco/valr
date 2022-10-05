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

template<class T, typename K>
void findClosestIvls(const IntervalTree<T, K>& tree, int start, int stop,
                 ivl_vector_t& closest,
                 std::pair<int, ivl_vector_t>& min_dist) {

  ivl_vector_t intervals = tree.intervals;
  typedef ivl_t interval;
  if (!intervals.empty() && !(stop < intervals.front().start)) {
    for (typename ivl_vector_t::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
      const interval& interval = *i;
      if (interval.stop > start && interval.start < stop) {
        // adjacent intervals are considered non-overlappping
        closest.push_back(interval);
      } else if (stop <= interval.start) {
        // find distance on left
        int ivl_dist_l = interval.start - stop ;
        // if smaller than best min dist found update pair with dist and intervals
        if (ivl_dist_l < min_dist.first) {
          min_dist.first = ivl_dist_l ;
          min_dist.second.clear() ;
          min_dist.second.push_back(interval) ;
        } else if (ivl_dist_l == min_dist.first) {
          // if same dist append intervals
          min_dist.second.push_back(interval) ;
        }
      } else if (start >= interval.stop) {
        // find distance on right
        int ivl_dist_r = start - interval.stop ;
        // if smaller than best min dist found update pair with dist and intervals
        if (ivl_dist_r < min_dist.first) {
          min_dist.first = ivl_dist_r ;
          min_dist.second.clear() ;
          min_dist.second.push_back(interval) ;
        } else if (ivl_dist_r == min_dist.first) {
          // if same dist append interval
          min_dist.second.push_back(interval) ;
        }
      }
    }
  }  else if (!intervals.empty()  && (stop <= intervals.front().start)) {
    for (typename ivl_vector_t::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
      const interval& interval = *i;
      if (interval.start > intervals.front().start) {
        continue ;
      } else {
        // find distance on left
        int ivl_dist_l = interval.start - stop ;
        // if smaller than best min dist found update pair with dist and intervals
        if (ivl_dist_l <= min_dist.first) {
          min_dist.first = ivl_dist_l;
          min_dist.second.clear() ;
          min_dist.second.push_back(interval) ;
        } else if (ivl_dist_l == min_dist.first) {
          // if same dist append intervals
          min_dist.second.push_back(interval) ;
        }
      }
    }
  }


  if (tree.left && start <= tree.center) {
    findClosestIvls(*tree.left, start, stop, closest, min_dist);
  }

  if (tree.right && stop >= tree.center) {
    findClosestIvls(*tree.right, start, stop,closest, min_dist);
  }

  // Finally report all of the non-overlapping closest intervals, only if at a left_node
  if (!(tree.right && tree.left)) {
    closest.insert(closest.end(), min_dist.second.begin(), min_dist.second.end())  ;
  }
}

#endif
