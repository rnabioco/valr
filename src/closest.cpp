// closest.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

// get overlap distance
inline int ivl_overlap(int xs, int xe, int ys, int ye) {
  return std::min(xe, ye) - std::max(xs, ys);
}

// class to store run-length encoding
// l = lengths
// v = values
// s = start positions for each value
template <typename T>
class RLE {
 public:
  std::vector<int> l;
  std::vector<int> v;
  std::vector<int> s;
  RLE(const T& x) {
    // based on https://rosettacode.org/wiki/Run-length_encoding#C++
    auto first = x.begin();
    auto last = x.end();
    auto idx = std::size_t{0};

    while (first != last) {
      auto const value = *first++;
      auto count = std::size_t{1};
      while (first != last && *first == value) {
        ++count;
        ++first;
      }
      v.push_back(value);
      l.push_back(count);
      s.push_back(idx);
      idx += count;
    }
  }
};

// obtain indexes of sorted vector
std::vector<int> sort_indexes(const std::vector<int>& v) {
  size_t n = v.size();
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  // stable sort is more efficient here when vector has many repeated values
  std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
}

// binary search
std::vector<int> findClosestPos(const std::vector<int>& x, const std::vector<int>& breaks) {
  std::size_t n = x.size();
  std::vector<int> out(n);
  auto breaks_it = breaks.begin();
  auto breaks_end = breaks.end();
  for (std::size_t i = 0; i < n; ++i) {
    auto ubound = std::upper_bound(breaks_it, breaks_end, x[i]);
    out[i] = std::distance(breaks_it, ubound);
  }
  return out;
}

void findClosestIvls(const std::vector<int>& q_starts, const std::vector<int>& q_ends,
                     const std::vector<int>& s_starts, const std::vector<int>& s_ends,
                     const std::vector<int>& gi_x, const std::vector<int>& gi_y,
                     std::vector<int>& indices_x, std::vector<int>& indices_y,
                     std::vector<int>& distance_sizes) {
  auto nx = q_starts.size();
  auto ny = s_starts.size();

  // sort ends if needed, storing sort order if unsorted
  // starts are always sorted
  std::vector<int> s_end_srted(ny);
  std::vector<int> s_ord_ends;

  bool e_is_sorted = std::is_sorted(s_ends.begin(), s_ends.end());
  if (!e_is_sorted) {
    s_ord_ends = sort_indexes(s_ends);
    for (size_t j = 0; j < ny; ++j) {
      s_end_srted[j] = s_ends[s_ord_ends[j]];
    }
  } else {
    for (size_t j = 0; j < ny; ++j) {
      s_end_srted[j] = s_ends[j];
    }
  }

  // rle encode to keep matches
  RLE<std::vector<int>> sstr_rle(s_starts);
  RLE<std::vector<int>> send_rle(s_end_srted);
  int sstr_n = sstr_rle.v.size();

  std::vector<int> before(nx);
  std::vector<int> after(nx);

  // compare x-ivl ends to y-ivl starts
  // q_ends - 1
  std::vector<int> q_ends_minus1(nx);
  for (size_t i = 0; i < nx; ++i) {
    q_ends_minus1[i] = q_ends[i] - 1;
  }
  before = findClosestPos(q_ends_minus1, sstr_rle.v);
  // compare x-ivl starts to y-ivl ends
  after = findClosestPos(q_starts, send_rle.v);

  int left_dist, right_dist, rle_idx, row_idx, run_length, overlap, dist;
  bool after_is_na, before_is_na, min_is_before;

  // compare before and after to determine closest
  for (size_t j = 0; j < nx; j++) {
    if (before[j] >= sstr_n) {
      before[j] = NA_INTEGER;
    }

    if (after[j] == 0) {
      after[j] = NA_INTEGER;
    } else {
      after[j] -= 1;
    }

    after_is_na = before_is_na = false;
    if (before[j] == NA_INTEGER) {
      before_is_na = true;
      left_dist = NA_INTEGER;
    } else {
      left_dist = sstr_rle.v[before[j]] - q_ends[j];
    }

    if (after[j] == NA_INTEGER) {
      after_is_na = true;
      right_dist = NA_INTEGER;
    } else {
      right_dist = q_starts[j] - send_rle.v[after[j]];
    }

    if (after_is_na && before_is_na) {
      // edge cases with overlapping or adjacent ivls
      // handled with a separate bed_intersect call
      continue;
    } else if (after_is_na || before_is_na) {
      min_is_before = after_is_na;
    } else {
      min_is_before = left_dist < right_dist;
    }

    if (min_is_before) {
      // before is closest
      rle_idx = before[j];
      row_idx = sstr_rle.s[rle_idx];
      run_length = sstr_rle.l[rle_idx];
      overlap = ivl_overlap(q_starts[j], q_ends[j], s_starts[row_idx], s_ends[row_idx]);
      if (overlap < 0) {
        overlap -= 1;
        dist = s_ends[row_idx] <= q_starts[j] ? overlap : -overlap;
        for (int k = 0; k < run_length; ++k) {
          indices_x.push_back(gi_x[j]);
          indices_y.push_back(gi_y[row_idx]);
          distance_sizes.push_back(dist);
          ++row_idx;
        }
      }

    } else {
      // handle same dist in before and after ivls
      if (left_dist == right_dist) {
        rle_idx = before[j];
        row_idx = sstr_rle.s[rle_idx];
        run_length = sstr_rle.l[rle_idx];
        overlap = ivl_overlap(q_starts[j], q_ends[j], s_starts[row_idx], s_ends[row_idx]);
        if (overlap < 0) {
          overlap -= 1;
          dist = s_ends[row_idx] <= q_starts[j] ? overlap : -overlap;
          for (int k = 0; k < run_length; ++k) {
            indices_x.push_back(gi_x[j]);
            indices_y.push_back(gi_y[row_idx]);
            distance_sizes.push_back(dist);
            ++row_idx;
          }
        }
      }
      // after is closest
      rle_idx = after[j];
      row_idx = send_rle.s[rle_idx];
      run_length = send_rle.l[rle_idx];
      row_idx = e_is_sorted ? row_idx : s_ord_ends[row_idx];
      overlap = ivl_overlap(q_starts[j], q_ends[j], s_starts[row_idx], s_ends[row_idx]);
      if (overlap < 0) {
        overlap -= 1;
        dist = s_ends[row_idx] <= q_starts[j] ? overlap : -overlap;
        for (int k = 0; k < run_length; ++k) {
          indices_x.push_back(gi_x[j]);
          indices_y.push_back(gi_y[row_idx]);
          distance_sizes.push_back(dist);
          ++row_idx;
        }
      }
    }
  }
}

[[cpp11::register]]
cpp11::writable::data_frame closest_impl(cpp11::data_frame gdf_x, cpp11::data_frame gdf_y,
                                         cpp11::integers grp_idx_x, cpp11::integers grp_idx_y,
                                         std::string suffix_x, std::string suffix_y) {
  valr::GroupedDataFrame grouped_x(gdf_x);
  valr::GroupedDataFrame grouped_y(gdf_y);

  cpp11::data_frame df_x = grouped_x.data();
  cpp11::data_frame df_y = grouped_y.data();

  // for subsetting / return df
  std::vector<int> indices_x;
  std::vector<int> indices_y;
  std::vector<int> distance_sizes;

  // Pre-allocate for typical case (at least one result per x interval)
  int estimated_size = df_x.nrow();
  indices_x.reserve(estimated_size);
  indices_y.reserve(estimated_size);
  distance_sizes.reserve(estimated_size);

  int ng_x = grp_idx_x.size();
  int ng_y = grp_idx_y.size();

  if (ng_x != ng_y) {
    cpp11::stop("incompatible groups found between x and y dataframes");
  }

  // access the group .rows list
  cpp11::list grp_indices_x = grouped_x.indices();
  cpp11::list grp_indices_y = grouped_y.indices();

  // Get start/end columns from dataframes
  cpp11::doubles df_x_starts = df_x["start"];
  cpp11::doubles df_x_ends = df_x["end"];
  cpp11::doubles df_y_starts = df_y["start"];
  cpp11::doubles df_y_ends = df_y["end"];

  // iterate through groups manually rather than using GroupApply
  // this keeps relevant data in vectors rather than
  // in vectors of interval type

  for (int i = 0; i < ng_x; i++) {
    // get next row index to subset from x and y groups
    // convert from R to C index
    int shared_x_index = grp_idx_x[i] - 1;
    int shared_y_index = grp_idx_y[i] - 1;

    // subset the group lists
    cpp11::integers gi_x_r(grp_indices_x[shared_x_index]);
    cpp11::integers gi_y_r(grp_indices_y[shared_y_index]);

    size_t nx = gi_x_r.size();
    size_t ny = gi_y_r.size();

    if (nx == 0 || ny == 0) {
      continue;
    }

    // Convert to 0-indexed std::vectors
    std::vector<int> gi_x(nx);
    std::vector<int> gi_y(ny);
    for (size_t j = 0; j < nx; ++j) {
      gi_x[j] = gi_x_r[j] - 1;
    }
    for (size_t j = 0; j < ny; ++j) {
      gi_y[j] = gi_y_r[j] - 1;
    }

    // Subset the coordinate vectors
    std::vector<int> q_starts(nx);
    std::vector<int> q_ends(nx);
    std::vector<int> s_starts(ny);
    std::vector<int> s_ends(ny);

    for (size_t j = 0; j < nx; ++j) {
      q_starts[j] = static_cast<int>(df_x_starts[gi_x[j]]);
      q_ends[j] = static_cast<int>(df_x_ends[gi_x[j]]);
    }
    for (size_t j = 0; j < ny; ++j) {
      s_starts[j] = static_cast<int>(df_y_starts[gi_y[j]]);
      s_ends[j] = static_cast<int>(df_y_ends[gi_y[j]]);
    }

    findClosestIvls(q_starts, q_ends, s_starts, s_ends, gi_x, gi_y, indices_x, indices_y,
                    distance_sizes);
  }

  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df_x, indices_x);
  cpp11::writable::data_frame subset_y = valr::subset_dataframe(df_y, indices_y);

  valr::DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, suffix_x, false);

  // y names, data
  out.add_df(subset_y, suffix_y, true);

  out.names.push_back(".dist");
  cpp11::writable::integers dist_vec(distance_sizes.size());
  for (size_t i = 0; i < distance_sizes.size(); ++i) {
    dist_vec[i] = distance_sizes[i];
  }
  out.data.push_back(dist_vec);

  int nrows = indices_x.size();
  return out.build(nrows);
}
