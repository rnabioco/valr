// closest.cpp
//
// Copyright (C) 2016 - 2023 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

// get overlap distance
inline int ivl_overlap(int xs, int xe, int ys, int ye) {
  return std::min(xe, ye) - std::max(xs, ys) ;
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
// adapted from user Lukasz Wiklendt https://stackoverflow.com/a/12399290/6276041
std::vector<int> sort_indexes(const IntegerVector& v) {
  size_t n = v.size();
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  // stable sort is more efficient here when vector has many repeated values
  std::stable_sort(idx.begin(), idx.end(),
                  [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// binary search
std::vector<int> findClosestPos(const IntegerVector& x,
                             const std::vector<int>& breaks) {
  std::size_t n = x.size();
  std::vector<int> out(n);
  auto breaks_it = breaks.begin();
  auto breaks_end = breaks.end();
  for(std::size_t i = 0; i < n; ++i) {
    auto ubound = std::upper_bound(breaks_it, breaks_end, x[i]);
    out[i] = std::distance(breaks_it, ubound);
  }
  return out;
}

void findClosestIvls(const IntegerVector& q_starts,
                    const IntegerVector& q_ends,
                    const IntegerVector& s_starts,
                    const IntegerVector& s_ends,
                    const IntegerVector& gi_x,
                    const IntegerVector& gi_y,
                    std::vector<int>& indices_x,
                    std::vector<int>& indices_y,
                    std::vector<int>& distance_sizes) {

  auto nx = q_starts.size();
  auto ny = s_starts.size();

  // sort ends if needed, storing sort order if unsorted
  // starts are always sorted
  std::vector<int> s_end_srted(ny);
  std::vector<int> s_ord_ends;

  bool e_is_sorted = std::is_sorted(s_ends.begin(), s_ends.end());
  if(!e_is_sorted){
    s_ord_ends = sort_indexes(s_ends);
    for(int j = 0; j < ny; ++j){
      s_end_srted[j] = s_ends[s_ord_ends[j]];
    }
  } else {
    for(int j = 0; j < ny; ++j){
      s_end_srted[j] = s_ends[j];
    }
  }

  // rle encode to keep matches
  RLE<IntegerVector> sstr_rle(s_starts);
  RLE<std::vector<int>> send_rle(s_end_srted);
  int sstr_n = sstr_rle.v.size();

  std::vector<int> before(nx);
  std::vector<int> after(nx);

  // compare x-ivl ends to y-ivl starts
  before = findClosestPos(q_ends - 1, sstr_rle.v);
  // compare x-ivl starts to y-ivl ends
  after = findClosestPos(q_starts, send_rle.v);


  int ld, rd, ni, idx, idx_l, ol, d;
  bool a_is_na, b_is_na, min_is_before;

  // compare before and after to determine closest
  for(int j = 0; j < nx; j++){
    if(before[j] >= sstr_n){
      before[j] = NA_INTEGER;
    }

    if(after[j] == 0){
      after[j] = NA_INTEGER;
    } else {
      after[j] -= 1;
    }

    a_is_na = b_is_na = false;
    if(IntegerVector::is_na(before[j])){
      b_is_na = true;
      ld = NA_INTEGER;
    } else {
      ld = sstr_rle.v[before[j]] - q_ends[j];
    }

    if(IntegerVector::is_na(after[j])){
      a_is_na = true;
      rd = NA_INTEGER;
    } else {
      rd = q_starts[j] - send_rle.v[after[j]];
    }

    if(a_is_na && b_is_na){
      // edge cases with overlapping or adjacent ivls
      // handled with a separate bed_intersect call
      continue;
    } else if (a_is_na || b_is_na) {
      min_is_before = a_is_na;
    } else {
      min_is_before = ld < rd;
    }

    if(min_is_before){
      // before is closest
      ni = before[j];
      idx = sstr_rle.s[ni];
      idx_l = sstr_rle.l[ni];
      ol = ivl_overlap(q_starts[j], q_ends[j],
                       s_starts[idx], s_ends[idx]);
      if(ol < 0) {
        ol -= 1;
        d = s_ends[idx] <= q_starts[j] ? ol : -ol;
        for (int k = 0; k < idx_l; ++k){
          indices_x.push_back(gi_x[j]);
          indices_y.push_back(gi_y[idx]);
          distance_sizes.push_back(d);
          ++idx;
        }
      }

    } else {
      // handle same dist in before and after ivls
      if(ld == rd){
        ni = before[j];
        idx = sstr_rle.s[ni];
        idx_l = sstr_rle.l[ni];
        ol = ivl_overlap(q_starts[j], q_ends[j],
                         s_starts[idx], s_ends[idx]);
        if(ol < 0) {
          ol -= 1;
          d = s_ends[idx] <= q_starts[j] ? ol : -ol;
          for (int k = 0; k < idx_l; ++k){
            indices_x.push_back(gi_x[j]);
            indices_y.push_back(gi_y[idx]);
            distance_sizes.push_back(d);
            ++idx;
          }
        }
      }
      // after is closest
      ni = after[j];
      idx = send_rle.s[ni];
      idx_l = send_rle.l[ni];
      idx = e_is_sorted ? idx : s_ord_ends[idx];
      ol = ivl_overlap(q_starts[j], q_ends[j], s_starts[idx], s_ends[idx]);
      if(ol < 0) {
        ol -= 1;
        d = s_ends[idx] <= q_starts[j] ? ol : -ol;
        for (int k = 0; k < idx_l; ++k){
          indices_x.push_back(gi_x[j]);
          indices_y.push_back(gi_y[idx]);
          distance_sizes.push_back(d);
          ++idx;
        }
      }
    }
  }

}


//[[Rcpp::export]]
DataFrame closest_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y,
                       IntegerVector grp_idx_x,
                       IntegerVector grp_idx_y,
                       const std::string& suffix_x,
                       const std::string& suffix_y) {

  DataFrame df_x = x.data() ;
  DataFrame df_y = y.data() ;

  // for subsetting / return df
  std::vector<int> indices_x ;
  std::vector<int> indices_y ;
  std::vector<int> distance_sizes ;

  int ng_x = grp_idx_x.size() ;
  int ng_y = grp_idx_y.size() ;

  if (ng_x != ng_y) {
    stop("incompatible groups found between x and y dataframes") ;
  }

  // access the group .rows list
  ListView grp_indices_x(x.indices()) ;
  ListView grp_indices_y(y.indices()) ;

  // iterate through groups manually rather than using GroupApply
  // this keeps relevant data in vectors rather than
  // in vectors of interval type

  for (int i = 0; i < ng_x; i++) {
    // get next row index to subset from x and y groups
    // convert from R to C index
    int shared_x_index = grp_idx_x[i] - 1;
    int shared_y_index = grp_idx_y[i] - 1;

    // subset the group lists
    IntegerVector gi_x, gi_y ;
    gi_x = grp_indices_x[shared_x_index];
    gi_x = gi_x - 1;
    gi_y = grp_indices_y[shared_y_index];
    gi_y = gi_y - 1;

    size_t nx = gi_x.size();
    size_t ny = gi_y.size();

    if (nx == 0 || ny == 0) {
      continue ;
    }

    IntegerVector q_starts = df_x["start"] ;
    IntegerVector q_ends   = df_x["end"] ;
    IntegerVector s_starts = df_y["start"];
    IntegerVector s_ends   = df_y["end"];

    q_starts = q_starts[gi_x];
    q_ends   = q_ends[gi_x];
    s_starts = s_starts[gi_y] ;
    s_ends   = s_ends[gi_y];

    findClosestIvls(q_starts, q_ends, s_starts, s_ends,
                    gi_x, gi_y,
                    indices_x, indices_y, distance_sizes);
  }

  DataFrame subset_x = subset_dataframe(df_x, indices_x) ;
  DataFrame subset_y = subset_dataframe(df_y, indices_y) ;

  DataFrameBuilder out;
  // x names, data
  out.add_df(subset_x, suffix_x, false) ;

  // y names, data
  out.add_df(subset_y, suffix_y, true) ;

  out.names.push_back(".dist") ;
  out.data.push_back(distance_sizes) ;

  auto nrows = subset_x.nrows() ;
  auto res = out.format_df(nrows) ;
  return res ;

}

