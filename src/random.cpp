// random.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

// [[Rcpp::export]]
DataFrame random_impl(DataFrame genome, int length, int n, int seed = 0) {

  CharacterVector chroms = genome["chrom"] ;
  NumericVector sizes = genome["size"] ;

  int nchrom = chroms.size() ;

  if (seed == 0)
    seed = round(R::runif(0, RAND_MAX)) ;

  // seed the generator
  auto generator = ENGINE(seed) ;

  // calculate weights for chrom distribution
  float mass = sum(sizes) ;
  NumericVector weights = sizes / mass ;

  Range chromidx(0, nchrom) ;
  PCONST_DIST chrom_dist(chromidx.begin(), chromidx.end(), weights.begin()) ;

  // make and store a DIST for each chrom size
  std::vector< UINT_DIST > size_rngs ;

  for (int i = 0; i < nchrom; ++i) {

    auto size = sizes[i] ;
    // sub length to avoid off-chrom coordinates
    UINT_DIST size_dist(0, size - length) ;
    size_rngs.push_back(size_dist) ;
  }

  CharacterVector rand_chroms(n) ;
  IntegerVector rand_starts(n) ;

  for (int i = 0; i < n; ++i) {

    auto chrom_idx = chrom_dist(generator) ;
    rand_chroms[i] = chroms[chrom_idx] ;

    UINT_DIST size_dist = size_rngs[chrom_idx] ;

    auto rand_start = size_dist(generator) ;
    rand_starts[i] = rand_start ;
  }

  IntegerVector rand_ends = rand_starts + length ;

  return DataFrame::create(_("chrom") = rand_chroms,
                           _("start") = rand_starts,
                           _("end") = rand_ends,
                           _("stringsAsFactors") = false) ;

}

/***R
library(dplyr)
genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 191822,
  "chr2", 17127713,
  "chr3", 11923987
)

# show chrom disribution
random_impl(genome, length = 1000, n = 1e6, seed = 0) %>%
  group_by(chrom) %>% summarize(n = n())

library(microbenchmark)
microbenchmark(
  random_impl(genome, length = 1000, n = 1e6, seed = 0),
  times = 10
)
*/
