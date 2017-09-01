// utils.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

// [[Rcpp::export]]
std::vector<int> unique_ids_impl(DataFrame x) {

  std::vector<int> ids ;
  int group_id = 0 ;
  std::unordered_map<std::string, int> group_map ;

  CharacterVector x_chroms = x["chrom"] ;
  IntegerVector x_starts = x["start"] ;
  IntegerVector x_ends   = x["end"] ;

  for (int i=0; i < x.nrows(); i++) {

    auto chrom = as<std::string>(x_chroms[i]) ;
    int x_start = x_starts[i] ;
    int x_end  = x_ends[i] ;

    std::string map_id = chrom + std::to_string(x_start) + std::to_string(x_end) ;

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
