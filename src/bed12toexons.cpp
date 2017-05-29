// bed12toexons.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

std::vector<int> csv_values(std::string csv) {

  std::vector<int> values ;
  std::stringstream ss(csv) ;

  // https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
  while (ss.good()) {
    std::string substr ;
    getline(ss, substr, ',') ;
    values.push_back(atoi(substr.c_str())) ;
  }

  return values ;
}

// [[Rcpp::export]]
DataFrame bed12toexons_impl(DataFrame x) {

  // input
  IntegerVector starts = x["start"] ;
  std::vector<std::string> exon_sizes = x["exon_sizes"] ;
  std::vector<std::string> exon_starts = x["exon_starts"] ;
  std::vector<std::string> strands = x["strand"] ;

  // storage
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
  std::vector<int> nums_out ; // exon numbers
  std::vector<int> df_idx ;

  // for looping
  std::vector<int> starts_exon ;
  std::vector<int> ends_exon ;
  std::vector<int> nums_exon ;

  for (int i = 0; i < starts.size(); i++) {

    std::vector<int> exon_start = csv_values(exon_starts[i]) ;
    std::vector<int> exon_size = csv_values(exon_sizes[i]) ;

    // calculate starts and ends for each exon
    int start = starts[i] ;
    for (int j = 0; j < exon_start.size(); j++) {
      starts_exon.push_back(start + exon_start[j]) ;
      ends_exon.push_back(start + exon_start[j] + exon_size[j]) ;
      df_idx.push_back(i) ;
    }

    // define range of exon numbers and reverse for `-` strand
    for (int k = 1; k < exon_start.size(); k++) {
      nums_exon.push_back(k) ;
    }

    if (strands[i] == "-") {
      std::reverse(nums_exon.begin(), nums_exon.end()) ;
    }

    // accumulate values
    starts_out.insert(starts_out.end(), starts_exon.begin(), starts_exon.end()) ;
    ends_out.insert(ends_out.end(), ends_exon.begin(), ends_exon.end()) ;
    nums_out.insert(nums_out.end(), nums_exon.begin(), nums_exon.end()) ;

    // reset local values
    starts_exon.clear() ;
    ends_exon.clear() ;
    nums_exon.clear() ;
  }

  CharacterVector dfnames = CharacterVector::create("chrom", "start", "end", "name", "score", "strand") ;
  DataFrame out = DataFrameSubsetVisitors(x, dfnames).subset(df_idx, "data.frame");

  out["start"] = starts_out ;
  out["end"] = ends_out ;
  out["score"] = nums_out ;

  return out ;
}

/***R
library(valr)
library(dplyr)
x <- read_bed12(valr_example('mm9.bed12.gz'))
devtools::load_all()
bed12toexons_impl(x) %>% as_data_frame()
*/
