// shuffle.cpp
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

typedef std::unordered_map<std::string, ivl_tree_t> chrom_tree_t ;
typedef std::unordered_map<std::string, ivl_vector_t> interval_map_t ;
typedef std::unordered_map<std::string, PCONST_DIST > interval_rng_t ;
typedef std::unordered_map<std::string, std::vector< UINT_DIST> > start_rng_t ;

chrom_tree_t makeIntervalTrees(DataFrame incl, interval_map_t interval_map) {

  chrom_tree_t chrom_trees ;

  // now build a map of chrom to interval tree
  for (auto kv : interval_map) {

    std::string chrom = kv.first ;
    ivl_vector_t iv = kv.second ;

    chrom_trees[chrom] = ivl_tree_t(std::move(iv)) ;

  }
  return chrom_trees ;
}

// used to select a chrom by its weighted mass
PCONST_DIST makeChromRNG(DataFrame incl) {

  CharacterVector incl_chroms = incl["chrom"] ;
  IntegerVector incl_starts = incl["start"] ;
  IntegerVector incl_ends = incl["end"] ;

  IntegerVector incl_sizes = incl_ends - incl_starts ;

  // keys sorted in lexographic order
  std::map<std::string, float> chrom_mass ;

  int nr = incl.nrows() ;
  for (int i = 0; i < nr; i++) {

    std::string chrom = as<std::string>(incl_chroms[i]) ;

    if (!chrom_mass.count(chrom))
      chrom_mass[chrom] = 0 ;

    float curr_mass = incl_sizes[i] ;
    chrom_mass[chrom] += curr_mass ;
  }

  std::vector<float> weights ;
  float total_mass = 0 ;
  for (auto kv : chrom_mass)  {
    auto mass = kv.second ;
    weights.push_back(mass) ;
    total_mass += mass ;
  }

  std::transform(weights.begin(), weights.end(), weights.begin(),
  [total_mass](float mass) {
    return mass / total_mass;
  }) ;

  CharacterVector chrom_names = unique(incl_chroms) ;
  int nchrom = chrom_names.size() ;

  // if there is a single chrom, then we can only sample the range [0,0]
  if (nchrom == 1) nchrom = 0;

  Range chrom_range(0, nchrom) ;
  PCONST_DIST chrom_rng(chrom_range.begin(), chrom_range.end(), weights.begin()) ;

  return chrom_rng ;
}

interval_map_t makeIntervalMap(DataFrame incl) {

  CharacterVector incl_chroms = incl["chrom"] ;
  IntegerVector incl_starts   = incl["start"] ;
  IntegerVector incl_ends     = incl["end"] ;

  int nr = incl.nrows() ;
  interval_map_t interval_map ;

  for (int i = 0; i < nr; ++i) {
    auto chrom = as<std::string>(incl_chroms[i]) ;

    if (!interval_map.count(chrom))
      interval_map[chrom] = ivl_vector_t() ;

    interval_map[chrom].push_back(ivl_t(incl_starts[i], incl_ends[i], i)) ;
  }

  return interval_map ;
}

// used to select an interval for a specific chrom
interval_rng_t makeIntervalWeights(interval_map_t interval_map) {

  interval_rng_t interval_map_rngs ;

  for (auto kv : interval_map) {

    auto chrom = kv.first ;
    auto intervals = kv.second ;

    float total_mass = 0 ;
    std::vector<float> weights ;

    for (auto i : intervals) {
      float mass = i.stop - i.start ;
      weights.push_back(mass) ;
      total_mass += mass ;
    }

    std::transform(weights.begin(), weights.end(), weights.begin(),
    [total_mass](float mass) {
      return mass / total_mass;
    }) ;

    auto n_ivls = intervals.size() ;

    // if there is a single interval, then we can only sample the range [0,0]
    if (n_ivls == 1) n_ivls = 0;

    Range ivl_range(0, n_ivls) ;
    PCONST_DIST ivl_rng(ivl_range.begin(), ivl_range.end(), weights.begin()) ;

    interval_map_rngs[chrom] = ivl_rng ;
  }

  return interval_map_rngs ;
}


start_rng_t makeStartRNGs(interval_map_t interval_map) {

  start_rng_t start_rngs ;

  for (auto kv : interval_map) {

    auto chrom = kv.first ;
    auto intervals = kv.second ;

    if (!start_rngs.count(chrom)) start_rngs[chrom] = { };

    for (auto i : intervals) {
      UINT_DIST rng(i.start, i.stop) ;
      start_rngs[chrom].push_back(rng) ;
    }
  }

  return start_rngs ;
}

//[[Rcpp::export]]
DataFrame shuffle_impl(DataFrame df, DataFrame incl, bool within = false,
                       int max_tries = 1000, int seed = 0) {

  // seed for reproducible intervals
  if (seed == 0) seed = round(R::runif(0, RAND_MAX)) ;

  // seed the generator
  auto generator = ENGINE(seed) ;

  // data on incoming df
  CharacterVector df_chroms = df["chrom"] ;
  IntegerVector df_starts   = df["start"] ;
  IntegerVector df_ends     = df["end"] ;

  IntegerVector df_sizes = df_ends - df_starts ;

  // RNG weighted by chromosome masses
  auto chrom_rng = makeChromRNG(incl) ;
  // map of chrom to intervals
  auto interval_map = makeIntervalMap(incl) ;
  // maps chroms to RNGs for interval index positions
  auto interval_rngs = makeIntervalWeights(interval_map) ;
  // maps chroms to RNGs for start dists
  auto start_rngs = makeStartRNGs(interval_map) ;
  // make a map of chrom to interval tree for each set of included intervals
  auto interval_trees = makeIntervalTrees(incl, interval_map) ;

  // storage for output
  int nr = df.nrows() ;
  CharacterVector chroms_out(nr) ;
  IntegerVector   starts_out(nr) ;
  IntegerVector   ends_out(nr) ;

  CharacterVector incl_chroms = incl["chrom"] ;
  // sort in lexographic order
  CharacterVector chrom_names = unique(incl_chroms).sort() ;

  for (int i = 0; i < nr; ++i) {

    // select a chromosome
    if (within) {
      chroms_out[i] = df_chroms[i] ;
    } else {
      // pick a random chrom index.
      int rand_idx = chrom_rng(generator) ;
      chroms_out[i] = chrom_names[rand_idx] ;
    }

    // get tree from map
    auto chrom = as<std::string>(chroms_out[i]) ;
    auto chrom_tree = interval_trees[chrom] ;

    bool inbounds = false ;
    int niter = 0 ;

    // get the interval rng
    auto interval_rng = interval_rngs[chrom] ;

    while (!inbounds) {

      niter++ ;
      if (niter > max_tries) {
        // tried too many times to find an overlap, bail
        stop("maximum iterations exceeded in bed_shuffle") ;
      }

      // get a random interval index
      int rand_ivl_idx = interval_rng(generator) ;
      // get the start rng and pick a start
      UINT_DIST start_rng = start_rngs[chrom][rand_ivl_idx] ;
      int rand_start = start_rng(generator) ;

      auto rand_end = rand_start + df_sizes[i] ;

      auto overlaps = chrom_tree.findOverlapping(rand_start, rand_end) ;

      // didn't find an overlap, keep going
      if (overlaps.empty()) continue ;

      // check that the chosen end is <= the end of the overlapping interval
      bool enclosed = true ;
      for (auto j : overlaps) {
        if (rand_start >= j.start) {
          if (rand_end > j.stop) {
            enclosed = false ;
          }
        }
      }
      if (!enclosed) continue ;

      // if we get here, all checks pass. keep the interval.
      inbounds = true ;

      starts_out[i] = rand_start ;
      ends_out[i] = rand_end ;
    }

  }

  return DataFrame::create(_("chrom") = chroms_out,
                           _("start") = starts_out,
                           _("end") = ends_out,
                           _("stringsAsFactors") = false) ;
}

/***R
library(dplyr)
library(valr)
library(testthat)
library(microbenchmark)

genome <- tibble::tribble(
   ~chrom,  ~size,
   "chr1", 50000000,
   "chr2", 60000000,
   "chr3", 80000000
)

incl <- tibble::tribble(
   ~chrom, ~start, ~end,
   "chr1", 1, 5000000,
   "chr1", 5000000, 50000000,
   "chr2", 1, 60000000,
   "chr3", 1, 80000000
)

x <- bed_random(genome, n = 100) %>% bed_sort()

shuffle_impl(x, incl) %>%
  group_by(chrom) %>%
  summarize(count = n())

library(microbenchmark)
# microbenchmark(shuffle_impl(x, incl), n = 10, unit = 's')

*/
