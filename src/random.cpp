// random.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

[[cpp11::register]]
writable::data_frame random_impl(data_frame genome, int length, int n, int seed = 0) {

  strings chroms = genome["chrom"] ;
  doubles sizes = genome["size"] ;

  int nchrom = chroms.size() ;

  if (seed == 0)
    seed = round(Rf_runif(0, RAND_MAX)) ;

  // seed the generator
  auto generator = ENGINE(seed) ;

  // calculate weights for chrom distribution
  double mass = std::accumulate(sizes.begin(), sizes.end(), 0.0); ;

  std::vector<double> weights(nchrom) ;
  for (int i = 0; i < nchrom; ++i) {
    weights[i] = sizes[i] / mass ;
  }

  std::vector<int> chromidx;
  for (int i = 0; i < nchrom; ++i) {
    chromidx.push_back(i);
  }

  PCONST_DIST chrom_dist(chromidx.begin(), chromidx.end(), weights.begin()) ;

  // make and store a DIST for each chrom size
  std::vector< UINT_DIST > size_rngs ;

  for (int i = 0; i < nchrom; ++i) {

    auto size = sizes[i] ;
    // sub length to avoid off-chrom coordinates
    UINT_DIST size_dist(0, size - length) ;
    size_rngs.push_back(size_dist) ;
  }

  std::vector<std::string> rand_chroms(n) ;
  std::vector<int> rand_starts(n) ;

  for (int i = 0; i < n; ++i) {

    int chrom_idx = chrom_dist(generator) ;
    rand_chroms[i] = chroms[chrom_idx] ;

    UINT_DIST size_dist = size_rngs[chrom_idx] ;

    auto rand_start = size_dist(generator) ;
    rand_starts[i] = rand_start ;
  }

  std::vector<int> rand_ends(rand_starts.size()) ;
  for (int i = 0; i < rand_starts.size(); ++i) {
    rand_ends[i] = rand_starts[i] + length ;
  }

  return writable::data_frame({
    "chrom"_nm = rand_chroms,
    "start"_nm = rand_starts,
    "end"_nm = rand_ends,
  });

}
