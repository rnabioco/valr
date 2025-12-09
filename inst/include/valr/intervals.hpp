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

// Simple sorted interval list with linear scan for overlaps
// Uses min_overlap parameter to control overlap semantics:
//   min_overlap = 1 (default): strict half-open, book-ended intervals do NOT overlap
//   min_overlap = 0: inclusive, book-ended intervals DO overlap (legacy behavior)
template <typename Coord = int, typename Value = int>
class IntervalTree {
 public:
  using interval = Interval<Coord, Value>;
  using interval_vector = std::vector<interval>;

 private:
  interval_vector intervals_;

 public:
  IntervalTree() {}

  // Build from interval vector (consumes the vector via move)
  explicit IntervalTree(interval_vector&& intervals) : intervals_(std::move(intervals)) {
    // Sort by start position for potential early termination
    std::sort(intervals_.begin(), intervals_.end(), IntervalStartCmp<Coord, Value>());
  }

  // Find all intervals overlapping [start, stop) with at least min_overlap bases
  // min_overlap = 1 (default): strict half-open (matches bedtools intersect, bedder-rs)
  // min_overlap = 0: book-ended intervals count as overlapping (legacy valr behavior)
  interval_vector findOverlapping(Coord start, Coord stop, Coord min_overlap = 1) const {
    interval_vector result;
    if (intervals_.empty())
      return result;

    result.reserve(16);  // Pre-allocate for typical case

    for (const auto& ivl : intervals_) {
      // Early termination depends on min_overlap:
      // - min_overlap > 0: can terminate when ivl.start >= stop (no overlap possible)
      // - min_overlap = 0: must include book-ended, so terminate when ivl.start > stop
      if (min_overlap > 0) {
        if (ivl.start >= stop)
          break;
      } else {
        if (ivl.start > stop)
          break;
      }

      // Calculate overlap: min(end1, end2) - max(start1, start2)
      Coord overlap = std::min(ivl.stop, stop) - std::max(ivl.start, start);
      if (overlap >= min_overlap) {
        result.push_back(ivl);
      }
    }
    return result;
  }

  // Visit all intervals overlapping [start, stop) with at least min_overlap bases
  template <typename Func>
  void visit_overlapping(Coord start, Coord stop, Func&& f, Coord min_overlap = 1) const {
    for (const auto& ivl : intervals_) {
      // Early termination depends on min_overlap
      if (min_overlap > 0) {
        if (ivl.start >= stop)
          break;
      } else {
        if (ivl.start > stop)
          break;
      }

      Coord overlap = std::min(ivl.stop, stop) - std::max(ivl.start, start);
      if (overlap >= min_overlap) {
        f(ivl);
      }
    }
  }

  bool empty() const { return intervals_.empty(); }
  size_t size() const { return intervals_.size(); }

  // Visit all intervals
  template <typename Func>
  void visit_all(Func&& f) const {
    for (const auto& ivl : intervals_) {
      f(ivl);
    }
  }
};

// Type aliases for common usage
using ivl_t = Interval<int, int>;
using ivl_vector_t = std::vector<ivl_t>;
using ivl_tree_t = IntervalTree<int, int>;

// Helper to build intervals from cpp11 vectors (uses doubles for coordinates)
inline ivl_vector_t make_intervals(const cpp11::doubles& starts, const cpp11::doubles& ends,
                                   const cpp11::integers& indices) {
  ivl_vector_t result;
  int n = starts.size();
  result.reserve(n);

  for (int i = 0; i < n; ++i) {
    result.emplace_back(static_cast<int>(starts[i]), static_cast<int>(ends[i]), indices[i]);
  }
  return result;
}

// Build interval vector from dataframe and group indices (1-based from R)
inline ivl_vector_t makeIntervalVector(const cpp11::data_frame& df,
                                       const cpp11::integers& indices) {
  cpp11::doubles starts = df["start"];
  cpp11::doubles ends = df["end"];

  ivl_vector_t ivls;
  int n = indices.size();
  ivls.reserve(n);

  for (int i = 0; i < n; ++i) {
    int j = indices[i] - 1;  // Convert from 1-based to 0-based
    ivls.emplace_back(static_cast<int>(starts[j]), static_cast<int>(ends[j]), j);
  }
  return ivls;
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
