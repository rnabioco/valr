/* Copyright (c) 2011 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

// Modified from IntervalTree.h from EKG
// https://github.com/ekg/intervaltree

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

template <class T, typename K = int>
class Interval {
public:
  K start;
  K stop;
  T value;
  Interval(K s, K e, const T& v)
    : start(s)
    , stop(e)
    , value(v)
  { }
};

template <class T, typename K>
K intervalStart(const Interval<T, K>& i) {
  return i.start;
}

template <class T, typename K>
K intervalStop(const Interval<T, K>& i) {
  return i.stop;
}

template <class T, typename K>
std::ostream& operator<<(std::ostream& out, Interval<T, K>& i) {
  out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
  return out;
}


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

template <class T, typename K = int>
class IntervalTree {

public:
  typedef Interval<T, K> interval;
  typedef std::vector<interval> intervalVector;
  typedef IntervalTree<T, K> intervalTree;

  intervalVector intervals;
  std::unique_ptr<intervalTree> left;
  std::unique_ptr<intervalTree> right;
  K center;

  IntervalTree<T, K>(void)
    : left(nullptr)
    , right(nullptr)
    , center(0)
  { }

private:
  std::unique_ptr<intervalTree> copyTree(const intervalTree& orig) {
    return std::unique_ptr<intervalTree>(new intervalTree(orig));
  }

public:

  IntervalTree<T, K>(const intervalTree& other)
    :   intervals(other.intervals),
        left(other.left ? copyTree(*other.left) : nullptr),
        right(other.right ? copyTree(*other.right) : nullptr),
        center(other.center)
  {
  }

public:

  IntervalTree<T, K>& operator=(const intervalTree& other) {
    center = other.center ;
    intervals = other.intervals;
    left = other.left ? copyTree(*other.left) : nullptr;
    right = other.right ? copyTree(*other.right) : nullptr;
    return *this;
  }

  // Note: changes the order of ivals
  IntervalTree<T, K>(
    intervalVector& ivals,
    std::size_t depth = 16,
    std::size_t minbucket = 64,
    K leftextent = 0,
    K rightextent = 0,
    std::size_t maxbucket = 512
  )
    : left(nullptr)
    , right(nullptr)
  {

    --depth;
    IntervalStartSorter<T, K> intervalStartSorter;
    if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
      std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
      intervals = ivals;
      center = ivals.at(ivals.size() / 2).start;
    } else {
      if (leftextent == 0 && rightextent == 0) {
        // sort intervals by start
        std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
      }

      K leftp = 0;
      K rightp = 0;
      K centerp = 0;

      if (leftextent || rightextent) {
        leftp = leftextent;
        rightp = rightextent;
      } else {
        leftp = ivals.front().start;
        std::vector<K> stops;
        stops.resize(ivals.size());
        transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T, K>);
        rightp = *max_element(stops.begin(), stops.end());
      }

      //centerp = ( leftp + rightp ) / 2;
      centerp = ivals.at(ivals.size() / 2).start;
      center = centerp;

      intervalVector lefts;
      intervalVector rights;

      for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i) {
        const interval& interval = *i;
        if (centerp > interval.stop) {
          lefts.push_back(interval);
        } else if (interval.start > centerp) {
          rights.push_back(interval);
        } else {
          intervals.push_back(interval);
        }
      }

      if (!lefts.empty()) {
        left = std::unique_ptr<intervalTree>(new intervalTree(lefts, depth, minbucket, leftp, centerp));
      }
      if (!rights.empty()) {
        right = std::unique_ptr<intervalTree>(new intervalTree(rights, depth, minbucket, centerp, rightp));
      }
    }
  }

  intervalVector findOverlapping(K start, K stop) const {
    intervalVector ov;
    this->findOverlapping(start, stop, ov);
    return ov;
  }

  void findOverlapping(K start, K stop, intervalVector& overlapping) const {
    if (!intervals.empty() && !(stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
        const interval& interval = *i;
        if (interval.stop >= start && interval.start <= stop) {
          overlapping.push_back(interval);
        }
      }
    }

    if (left && start <= center) {
      left->findOverlapping(start, stop, overlapping);
    }

    if (right && stop >= center) {
      right->findOverlapping(start, stop, overlapping);
    }

  }

  intervalVector findContained(K start, K stop) const {
    intervalVector contained;
    this->findContained(start, stop, contained);
    return contained;
  }

  void findContained(K start, K stop, intervalVector& contained) const {
    if (!intervals.empty() && !(stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
        const interval& interval = *i;
        if (interval.start >= start && interval.stop <= stop) {
          contained.push_back(interval);
        }
      }
    }

    if (left && start <= center) {
      left->findContained(start, stop, contained);
    }

    if (right && stop >= center) {
      right->findContained(start, stop, contained);
    }

  }

  intervalVector findClosest(K start, K stop) const {
    intervalVector closest ;
    std::pair<int, intervalVector> min_dist_l, min_dist_r;
    this->findClosest(start, stop, closest, min_dist_l, min_dist_r) ;
    return closest;
  }

  void findClosest(K start, K stop, intervalVector&  closest,
                   std::pair<int, intervalVector>& min_dist) const {

    if (!intervals.empty() && !(stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
        const interval& interval = *i;
        if (interval.stop >= start && interval.start <= stop) {
          closest.push_back(interval);
        } else if (stop < interval.start) {
          // finddistance on left
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
        } else if (start > interval.stop) {
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
    }  else if (!intervals.empty()  && (stop < intervals.front().start)) {
      for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
        const interval& interval = *i;
        if (interval.start > intervals.front().start) {
          continue ;
        } else {
          // finddistance on left
          int ivl_dist_l = interval.start - stop ;
          // if smaller than best min dist found update pair with dist and intervals
          if (ivl_dist_l < min_dist.first) {
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


    if (left && start <= center) {
      left->findClosest(start, stop,  closest, min_dist);
    }

    if (right && stop >= center) {
      right->findClosest(start, stop,  closest, min_dist);
    }

    // Finally report all of the non-overlapping closest intervals, only if at a left_node
    if (!(right && left)) {
      closest.insert(closest.end(), min_dist.second.begin(), min_dist.second.end())  ;
    }
  }
  ~IntervalTree(void) = default;

};

#endif
