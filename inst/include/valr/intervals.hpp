// intervals.hpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef VALR_INTERVALS_HPP
#define VALR_INTERVALS_HPP

#include <cpp11.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace valr {

// Basic interval structure
template <typename Coord = int, typename Value = int>
struct Interval {
  Coord start;
  Coord stop;  // Using 'stop' for compatibility with existing code
  Value value;

  Interval() : start(0), stop(0), value(0) {}
  Interval(Coord s, Coord e, Value v) : start(s), stop(e), value(v) {}
};

// Compute overlap between two intervals
template <typename Coord>
inline Coord interval_overlap(Coord s1, Coord e1, Coord s2, Coord e2) {
  return std::min(e1, e2) - std::max(s1, s2);
}

template <typename Coord, typename Value>
inline Coord interval_overlap(const Interval<Coord, Value>& a, const Interval<Coord, Value>& b) {
  return interval_overlap(a.start, a.stop, b.start, b.stop);
}

// Sort comparators
template <typename Coord, typename Value>
struct IntervalStartCmp {
  bool operator()(const Interval<Coord, Value>& a, const Interval<Coord, Value>& b) const {
    return a.start < b.start;
  }
};

template <typename Coord, typename Value>
struct IntervalStopCmp {
  bool operator()(const Interval<Coord, Value>& a, const Interval<Coord, Value>& b) const {
    return a.stop < b.stop;
  }
};

template <typename Coord, typename Value>
struct IntervalStopDescCmp {
  bool operator()(const Interval<Coord, Value>& a, const Interval<Coord, Value>& b) const {
    return a.stop > b.stop;
  }
};

// Optimized Interval Tree using implicit binary tree layout
// Uses augmented max-end for efficient pruning during queries
template <typename Coord = int, typename Value = int>
class IntervalTree {
 public:
  using interval = Interval<Coord, Value>;
  using interval_vector = std::vector<interval>;

 private:
  // Flat array storage for cache efficiency
  // Tree is stored in heap order: left(i) = 2*i+1, right(i) = 2*i+2
  std::vector<Coord> starts_;
  std::vector<Coord> stops_;
  std::vector<Value> values_;
  std::vector<Coord> max_stops_;  // Augmented: max stop in subtree rooted at i
  size_t size_;

  // Build augmented max values bottom-up
  void build_augmented() {
    if (size_ == 0)
      return;

    // Process nodes from bottom to top (reverse level order)
    for (int i = static_cast<int>(size_) - 1; i >= 0; --i) {
      Coord max_stop = stops_[i];

      size_t left = 2 * i + 1;
      size_t right = 2 * i + 2;

      if (left < size_) {
        max_stop = std::max(max_stop, max_stops_[left]);
      }
      if (right < size_) {
        max_stop = std::max(max_stop, max_stops_[right]);
      }

      max_stops_[i] = max_stop;
    }
  }

  // Recursive helper for building balanced tree
  void build_tree(interval_vector& intervals, size_t idx, size_t lo, size_t hi) {
    if (lo >= hi || idx >= size_)
      return;

    size_t mid = lo + (hi - lo) / 2;

    // Place median at current index
    starts_[idx] = intervals[mid].start;
    stops_[idx] = intervals[mid].stop;
    values_[idx] = intervals[mid].value;

    // Recursively build left and right subtrees
    if (mid > lo) {
      build_tree(intervals, 2 * idx + 1, lo, mid);
    }
    if (mid + 1 < hi) {
      build_tree(intervals, 2 * idx + 2, mid + 1, hi);
    }
  }

  // Query helper - recursively find overlapping intervals
  template <typename OutputIt>
  void query_overlapping(size_t idx, Coord qstart, Coord qstop, OutputIt& out) const {
    if (idx >= size_)
      return;

    // Prune: if max_stop in this subtree < query start, no overlaps possible
    if (max_stops_[idx] < qstart)
      return;

    // Check left subtree first (intervals with smaller starts)
    size_t left = 2 * idx + 1;
    if (left < size_) {
      query_overlapping(left, qstart, qstop, out);
    }

    // Check current interval
    if (starts_[idx] <= qstop && stops_[idx] >= qstart) {
      *out++ = interval(starts_[idx], stops_[idx], values_[idx]);
    }

    // Prune right subtree: if current start > query stop, no need to check right
    // (since all intervals in right subtree have start >= current start)
    if (starts_[idx] > qstop)
      return;

    // Check right subtree
    size_t right = 2 * idx + 2;
    if (right < size_) {
      query_overlapping(right, qstart, qstop, out);
    }
  }

 public:
  IntervalTree() : size_(0) {}

  // Build from interval vector (consumes the vector via move)
  explicit IntervalTree(interval_vector&& intervals) {
    size_ = intervals.size();
    if (size_ == 0)
      return;

    // Sort by start position
    std::sort(intervals.begin(), intervals.end(), IntervalStartCmp<Coord, Value>());

    // Allocate storage
    starts_.resize(size_);
    stops_.resize(size_);
    values_.resize(size_);
    max_stops_.resize(size_);

    // Build balanced tree
    build_tree(intervals, 0, 0, size_);

    // Build augmented max values
    build_augmented();
  }

  // Find all intervals overlapping [start, stop]
  interval_vector findOverlapping(Coord start, Coord stop) const {
    interval_vector result;
    if (size_ == 0)
      return result;

    result.reserve(16);  // Pre-allocate for typical case
    auto inserter = std::back_inserter(result);
    query_overlapping(0, start, stop, inserter);
    return result;
  }

  // Visit all intervals overlapping [start, stop]
  template <typename Func>
  void visit_overlapping(Coord start, Coord stop, Func&& f) const {
    if (size_ == 0)
      return;

    struct Visitor {
      Func& func;
      void operator=(const interval& ivl) { func(ivl); }
    };
    Visitor visitor{f};
    query_overlapping(0, start, stop, visitor);
  }

  bool empty() const { return size_ == 0; }
  size_t size() const { return size_; }

  // Visit all intervals in tree
  template <typename Func>
  void visit_all(Func&& f) const {
    for (size_t i = 0; i < size_; ++i) {
      f(interval(starts_[i], stops_[i], values_[i]));
    }
  }
};

// Type aliases for common usage
using ivl_t = Interval<int, int>;
using ivl_vector_t = std::vector<ivl_t>;
using ivl_tree_t = IntervalTree<int, int>;

// Helper to build intervals from cpp11 vectors
inline ivl_vector_t make_intervals(const cpp11::integers& starts, const cpp11::integers& ends,
                                   const cpp11::integers& indices) {
  ivl_vector_t result;
  int n = starts.size();
  result.reserve(n);

  for (int i = 0; i < n; ++i) {
    result.emplace_back(starts[i], ends[i], indices[i]);
  }
  return result;
}

// Legacy function name for compatibility
inline int intervalOverlap(const ivl_t& a, const ivl_t& b) {
  return interval_overlap(a, b);
}

// Legacy sorter for compatibility
template <typename Coord, typename Value>
using IntervalSorterDesc = IntervalStopDescCmp<Coord, Value>;

}  // namespace valr

// Global typedefs for backward compatibility during migration
using ivl_t = valr::ivl_t;
using ivl_vector_t = valr::ivl_vector_t;
using ivl_tree_t = valr::ivl_tree_t;

template <typename Coord, typename Value>
using IntervalSorterDesc = valr::IntervalStopDescCmp<Coord, Value>;

inline int intervalOverlap(const ivl_t& a, const ivl_t& b) {
  return valr::interval_overlap(a, b);
}

#endif  // VALR_INTERVALS_HPP
