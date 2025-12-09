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
#include <memory>
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

// Augmented interval tree using centered approach
// Inspired by Erik Garrison's implementation but simplified for valr's needs
//
// Key features:
// - O(log n + k) query time where k = number of overlaps
// - Supports min_overlap parameter for valr's overlap semantics
// - Falls back to linear scan for small datasets (< 64 intervals)
// - Memory efficient: reuses interval storage
template <typename Coord = int, typename Value = int>
class IntervalTree {
 public:
  using interval = Interval<Coord, Value>;
  using interval_vector = std::vector<interval>;

 private:
  // Intervals stored at this node (those spanning the center)
  interval_vector intervals_;
  // Subtrees
  std::unique_ptr<IntervalTree> left_;
  std::unique_ptr<IntervalTree> right_;
  // Center point for partitioning
  Coord center_;
  // Maximum endpoint in this subtree (for pruning)
  Coord max_stop_;

  // Threshold below which we use linear scan instead of tree.
  // For small interval sets, the overhead of tree construction and traversal
  // exceeds the benefit of O(log n) lookups. Empirically, 64 is a good cutoff.
  static constexpr size_t LINEAR_THRESHOLD = 64;

 public:
  IntervalTree() : center_(0), max_stop_(std::numeric_limits<Coord>::min()) {}

  // Build from interval vector (consumes the vector via move)
  explicit IntervalTree(interval_vector&& intervals, size_t depth = 16, size_t min_bucket = 64)
      : left_(nullptr), right_(nullptr), center_(0), max_stop_(std::numeric_limits<Coord>::min()) {
    if (intervals.empty()) {
      return;
    }

    // Find min/max for computing center
    Coord min_start = std::numeric_limits<Coord>::max();
    Coord max_stop = std::numeric_limits<Coord>::min();
    for (const auto& ivl : intervals) {
      min_start = std::min(min_start, ivl.start);
      max_stop = std::max(max_stop, ivl.stop);
    }
    max_stop_ = max_stop;

    // For small sets or at depth limit, store linearly
    if (depth == 0 || intervals.size() < min_bucket) {
      // Sort by start for efficient querying
      std::sort(intervals.begin(), intervals.end(), IntervalStartCmp<Coord, Value>());
      intervals_ = std::move(intervals);
      return;
    }

    // Choose center as midpoint of the range
    center_ = (min_start + max_stop) / 2;

    // Partition intervals into left, center, right
    interval_vector lefts, rights;
    lefts.reserve(intervals.size() / 2);
    rights.reserve(intervals.size() / 2);

    for (auto& ivl : intervals) {
      if (ivl.stop < center_) {
        // Entirely to the left of center
        lefts.push_back(std::move(ivl));
      } else if (ivl.start > center_) {
        // Entirely to the right of center
        rights.push_back(std::move(ivl));
      } else {
        // Spans the center - store at this node
        intervals_.push_back(std::move(ivl));
      }
    }

    // Sort center intervals by start for querying
    std::sort(intervals_.begin(), intervals_.end(), IntervalStartCmp<Coord, Value>());

    // Recursively build subtrees
    if (!lefts.empty()) {
      left_ = std::make_unique<IntervalTree>(std::move(lefts), depth - 1, min_bucket);
    }
    if (!rights.empty()) {
      right_ = std::make_unique<IntervalTree>(std::move(rights), depth - 1, min_bucket);
    }
  }

  // Find all intervals overlapping [start, stop) with at least min_overlap bases
  // min_overlap = 1 (default): strict half-open (matches bedtools)
  // min_overlap = 0: book-ended intervals count as overlapping (legacy valr)
  interval_vector findOverlapping(Coord start, Coord stop, Coord min_overlap = 1) const {
    interval_vector result;
    result.reserve(16);
    findOverlappingImpl(start, stop, min_overlap, result);
    return result;
  }

  // Visit all intervals overlapping [start, stop) with at least min_overlap bases
  template <typename Func>
  void visit_overlapping(Coord start, Coord stop, Func&& f, Coord min_overlap = 1) const {
    visitOverlappingImpl(start, stop, std::forward<Func>(f), min_overlap);
  }

  bool empty() const {
    if (!intervals_.empty())
      return false;
    if (left_ && !left_->empty())
      return false;
    if (right_ && !right_->empty())
      return false;
    return true;
  }

  size_t size() const {
    size_t count = intervals_.size();
    if (left_)
      count += left_->size();
    if (right_)
      count += right_->size();
    return count;
  }

  // Visit all intervals in the tree
  template <typename Func>
  void visit_all(Func&& f) const {
    if (left_)
      left_->visit_all(f);
    for (const auto& ivl : intervals_) {
      f(ivl);
    }
    if (right_)
      right_->visit_all(f);
  }

 private:
  void findOverlappingImpl(Coord start, Coord stop, Coord min_overlap,
                           interval_vector& result) const {
    // Check if we're a leaf node (no subtrees, just linear storage)
    if (!left_ && !right_) {
      // Linear scan through intervals_
      for (const auto& ivl : intervals_) {
        // Early termination: sorted by start, so if ivl.start is past our query, done
        if (min_overlap > 0) {
          if (ivl.start >= stop)
            break;
        } else {
          if (ivl.start > stop)
            break;
        }

        Coord overlap = std::min(ivl.stop, stop) - std::max(ivl.start, start);
        if (overlap >= min_overlap) {
          result.push_back(ivl);
        }
      }
      return;
    }

    // Check left subtree
    if (left_ && start < center_) {
      // Prune if query is entirely past left subtree's max endpoint
      if (min_overlap > 0) {
        if (start < left_->max_stop_) {
          left_->findOverlappingImpl(start, stop, min_overlap, result);
        }
      } else {
        if (start <= left_->max_stop_) {
          left_->findOverlappingImpl(start, stop, min_overlap, result);
        }
      }
    }

    // Check intervals at this node (they all span center_)
    for (const auto& ivl : intervals_) {
      Coord overlap = std::min(ivl.stop, stop) - std::max(ivl.start, start);
      if (overlap >= min_overlap) {
        result.push_back(ivl);
      }
    }

    // Check right subtree
    if (right_ && stop > center_) {
      right_->findOverlappingImpl(start, stop, min_overlap, result);
    }
  }

  template <typename Func>
  void visitOverlappingImpl(Coord start, Coord stop, Func&& f, Coord min_overlap) const {
    // Check if we're a leaf node (no subtrees, just linear storage)
    if (!left_ && !right_) {
      for (const auto& ivl : intervals_) {
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
      return;
    }

    // Check left subtree
    if (left_ && start < center_) {
      if (min_overlap > 0) {
        if (start < left_->max_stop_) {
          left_->visitOverlappingImpl(start, stop, f, min_overlap);
        }
      } else {
        if (start <= left_->max_stop_) {
          left_->visitOverlappingImpl(start, stop, f, min_overlap);
        }
      }
    }

    // Check intervals at this node
    for (const auto& ivl : intervals_) {
      Coord overlap = std::min(ivl.stop, stop) - std::max(ivl.start, start);
      if (overlap >= min_overlap) {
        f(ivl);
      }
    }

    // Check right subtree
    if (right_ && stop > center_) {
      right_->visitOverlappingImpl(start, stop, f, min_overlap);
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
