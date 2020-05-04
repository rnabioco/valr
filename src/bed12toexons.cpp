// bed12toexons.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
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

    if (substr.empty()) break ;

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

  for (int i = 0; i < starts.size(); i++) {

    std::vector<int> exon_start = csv_values(exon_starts[i]) ;
    std::vector<int> exon_size = csv_values(exon_sizes[i]) ;

    // calculate starts and ends for each exon
    int start = starts[i] ;
    int n = exon_start.size() ;

    for (int j = 0; j < n; j++) {

      starts_out.push_back(start + exon_start[j]) ;
      ends_out.push_back(start + exon_start[j] + exon_size[j]) ;

      if (strands[i] == "-") {
        nums_out.push_back(n - j) ;
      } else {
        nums_out.push_back(j + 1) ;
      }

      df_idx.push_back(i) ;
    }
  }

  DataFrame out = subset_dataframe(x, df_idx) ;

  out["start"] = starts_out ;
  out["end"] = ends_out ;
  out["score"] = nums_out ;

  return out ;
}

/***R
library(valr)
library(dplyr)
x <- read_bed12(valr_example('mm9.refGene.bed.gz'))
bed12_to_exons(x)
*/
