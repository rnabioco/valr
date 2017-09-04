// utils.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

// https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
inline void hash_combine(std::size_t& seed) { }

template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  hash_combine(seed, rest...);
}

// [[Rcpp::export]]
std::vector<int> unique_ids_impl(DataFrame x) {

  std::vector<int> ids ;
  int group_id = 0 ;
  std::unordered_map<size_t, int> group_map ;

  CharacterVector x_chroms = x["chrom"] ;
  IntegerVector x_starts = x["start"] ;
  IntegerVector x_ends   = x["end"] ;

  for (int i=0; i < x.nrows(); i++) {

    auto chrom = as<std::string>(x_chroms[i]) ;
    int x_start = x_starts[i] ;
    int x_end  = x_ends[i] ;

    size_t map_id = 0 ;
    hash_combine(map_id, chrom, x_start, x_end) ;

    auto it = group_map.find(map_id) ;
    if (it == group_map.end()) {
      // map_id has not been seen yet. increment and insert.
      group_id++ ;
      group_map.insert({map_id, group_id}) ;
    } else {
      // already seen so use that id
      group_id = it->second ;
    }

    ids.push_back(group_id) ;
  }

  return ids ;
}

// [[Rcpp::export]]
DataFrame unmatched_groups_impl(GroupedDataFrame x, GroupedDataFrame y) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;

  DataFrame labels_x(data_x.attr("labels"));
  DataFrame labels_y(data_y.attr("labels"));

  std::vector<int> indices ;

  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for (int nx = 0; nx < ng_x; nx++, ++git_x) {

    GroupedSlicingIndex gi_x = *git_x ;

    bool match = true ;
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;

    for (int ny = 0; ny < ng_y; ny++, ++git_y) {
      match = compare_rows(labels_x, labels_y, nx, ny);
      if (match) break ;
    }

    if (match) continue ;

    for (int i = 0; i < gi_x.size(); i++) {
      indices.push_back(gi_x[i]) ;
    }
  }

  DataFrame subset_x = DataFrameSubsetVisitors(data_x, data_x.names()).subset(indices, "data.frame");

  return(subset_x) ;

}
