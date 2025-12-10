// shuffle.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <Rmath.h>

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"
#include "valr/random.hpp"

using namespace cpp11::literals;

typedef std::unordered_map<std::string, valr::ivl_tree_t> chrom_tree_t;
typedef std::unordered_map<std::string, valr::ivl_vector_t> interval_map_t;
typedef std::unordered_map<std::string, valr::PiecewiseConstDist> interval_rng_t;
typedef std::unordered_map<std::string, std::vector<valr::UniformIntDist>> start_rng_t;

chrom_tree_t makeIntervalTrees(interval_map_t interval_map) {
  chrom_tree_t chrom_trees;

  // build a map of chrom to interval tree
  for (auto& kv : interval_map) {
    std::string chrom = kv.first;
    valr::ivl_vector_t iv = kv.second;
    chrom_trees[chrom] = valr::ivl_tree_t(std::move(iv));
  }
  return chrom_trees;
}

// used to select a chrom by its weighted mass
valr::PiecewiseConstDist makeChromRNG(cpp11::data_frame incl) {
  cpp11::strings incl_chroms = incl["chrom"];
  cpp11::doubles incl_starts = incl["start"];
  cpp11::doubles incl_ends = incl["end"];

  // keys sorted in lexographic order
  std::map<std::string, float> chrom_mass;

  int nr = incl.nrow();
  for (int i = 0; i < nr; i++) {
    std::string chrom(incl_chroms[i]);

    if (!chrom_mass.count(chrom))
      chrom_mass[chrom] = 0;

    float curr_mass = static_cast<float>(incl_ends[i] - incl_starts[i]);
    chrom_mass[chrom] += curr_mass;
  }

  std::vector<float> weights;
  float total_mass = 0;
  for (auto& kv : chrom_mass) {
    auto mass = kv.second;
    weights.push_back(mass);
    total_mass += mass;
  }

  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [total_mass](float mass) { return mass / total_mass; });

  // Count unique chroms
  std::set<std::string> unique_chroms;
  for (int i = 0; i < nr; i++) {
    unique_chroms.insert(std::string(incl_chroms[i]));
  }
  int nchrom = unique_chroms.size();

  // if there is a single chrom, then we can only sample the range [0,0]
  if (nchrom == 1)
    nchrom = 0;

  std::vector<int> chrom_range;
  for (int i = 0; i <= nchrom; i++) {
    chrom_range.push_back(i);
  }
  valr::PiecewiseConstDist chrom_rng(chrom_range.begin(), chrom_range.end(), weights.begin());

  return chrom_rng;
}

interval_map_t makeIntervalMap(cpp11::data_frame incl) {
  cpp11::strings incl_chroms = incl["chrom"];
  cpp11::doubles incl_starts = incl["start"];
  cpp11::doubles incl_ends = incl["end"];

  int nr = incl.nrow();
  interval_map_t interval_map;

  for (int i = 0; i < nr; ++i) {
    std::string chrom(incl_chroms[i]);

    if (!interval_map.count(chrom))
      interval_map[chrom] = valr::ivl_vector_t();

    int start = static_cast<int>(incl_starts[i]);
    int end = static_cast<int>(incl_ends[i]);
    interval_map[chrom].push_back(valr::ivl_t(start, end, i));
  }

  return interval_map;
}

// used to select an interval for a specific chrom
interval_rng_t makeIntervalWeights(interval_map_t interval_map) {
  interval_rng_t interval_map_rngs;

  for (auto& kv : interval_map) {
    auto chrom = kv.first;
    auto intervals = kv.second;

    float total_mass = 0;
    std::vector<float> weights;

    for (auto& i : intervals) {
      float mass = static_cast<float>(i.stop - i.start);
      weights.push_back(mass);
      total_mass += mass;
    }

    std::transform(weights.begin(), weights.end(), weights.begin(),
                   [total_mass](float mass) { return mass / total_mass; });

    auto n_ivls = intervals.size();

    // if there is a single interval, then we can only sample the range [0,0]
    if (n_ivls == 1)
      n_ivls = 0;

    std::vector<int> ivl_range;
    for (size_t i = 0; i <= n_ivls; i++) {
      ivl_range.push_back(static_cast<int>(i));
    }
    valr::PiecewiseConstDist ivl_rng(ivl_range.begin(), ivl_range.end(), weights.begin());

    interval_map_rngs[chrom] = ivl_rng;
  }

  return interval_map_rngs;
}

start_rng_t makeStartRNGs(interval_map_t interval_map) {
  start_rng_t start_rngs;

  for (auto& kv : interval_map) {
    auto chrom = kv.first;
    auto intervals = kv.second;

    if (!start_rngs.count(chrom))
      start_rngs[chrom] = {};

    for (auto& i : intervals) {
      valr::UniformIntDist rng(i.start, i.stop);
      start_rngs[chrom].push_back(rng);
    }
  }

  return start_rngs;
}

[[cpp11::register]] cpp11::writable::data_frame shuffle_impl(cpp11::data_frame df,
                                                             cpp11::data_frame incl,
                                                             bool within = false,
                                                             int max_tries = 1000, int seed = 0) {
  // seed for reproducible intervals
  if (seed == 0)
    seed = static_cast<int>(Rf_runif(0, RAND_MAX));

  // seed the generator
  auto generator = valr::Engine(seed);

  // data on incoming df
  cpp11::strings df_chroms = df["chrom"];
  cpp11::doubles df_starts = df["start"];
  cpp11::doubles df_ends = df["end"];

  int nr = df.nrow();
  std::vector<int> df_sizes(nr);
  for (int i = 0; i < nr; i++) {
    df_sizes[i] = static_cast<int>(df_ends[i] - df_starts[i]);
  }

  // RNG weighted by chromosome masses
  auto chrom_rng = makeChromRNG(incl);
  // map of chrom to intervals
  auto interval_map = makeIntervalMap(incl);
  // maps chroms to RNGs for interval index positions
  auto interval_rngs = makeIntervalWeights(interval_map);
  // maps chroms to RNGs for start dists
  auto start_rngs = makeStartRNGs(interval_map);
  // make a map of chrom to interval tree for each set of included intervals
  auto interval_trees = makeIntervalTrees(interval_map);

  // storage for output
  cpp11::writable::strings chroms_out(nr);
  cpp11::writable::integers starts_out(nr);
  cpp11::writable::integers ends_out(nr);

  cpp11::strings incl_chroms = incl["chrom"];

  // Get unique chroms sorted lexographically
  std::set<std::string> unique_chrom_set;
  for (int i = 0; i < incl.nrow(); i++) {
    unique_chrom_set.insert(std::string(incl_chroms[i]));
  }
  std::vector<std::string> chrom_names(unique_chrom_set.begin(), unique_chrom_set.end());

  for (int i = 0; i < nr; ++i) {
    std::string selected_chrom;

    // select a chromosome
    if (within) {
      selected_chrom = std::string(df_chroms[i]);
    } else {
      // pick a random chrom index.
      int rand_idx = chrom_rng(generator);
      selected_chrom = chrom_names[rand_idx];
    }
    chroms_out[i] = selected_chrom;

    // get tree from map - verify chrom exists
    auto tree_it = interval_trees.find(selected_chrom);
    if (tree_it == interval_trees.end()) {
      cpp11::stop("chromosome not found in inclusion regions: %s", selected_chrom.c_str());
    }
    auto& chrom_tree = tree_it->second;

    bool inbounds = false;
    int niter = 0;

    // get the interval rng - verify chrom exists
    auto rng_it = interval_rngs.find(selected_chrom);
    if (rng_it == interval_rngs.end()) {
      cpp11::stop("chromosome not found in interval RNGs: %s", selected_chrom.c_str());
    }
    auto interval_rng = rng_it->second;

    // get start RNGs for this chrom - verify exists
    auto start_rng_it = start_rngs.find(selected_chrom);
    if (start_rng_it == start_rngs.end()) {
      cpp11::stop("chromosome not found in start RNGs: %s", selected_chrom.c_str());
    }
    auto& chrom_start_rngs = start_rng_it->second;

    while (!inbounds) {
      niter++;
      if (niter > max_tries) {
        // tried too many times to find an overlap, bail
        cpp11::stop("maximum iterations exceeded in bed_shuffle");
      }

      // get a random interval index
      int rand_ivl_idx = interval_rng(generator);

      // bounds check on random interval index
      if (rand_ivl_idx < 0 || static_cast<size_t>(rand_ivl_idx) >= chrom_start_rngs.size()) {
        continue;  // invalid index, try again
      }

      // get the start rng and pick a start
      valr::UniformIntDist start_rng = chrom_start_rngs[rand_ivl_idx];
      int rand_start = start_rng(generator);

      int rand_end = rand_start + df_sizes[i];

      auto overlaps = chrom_tree.findOverlapping(rand_start, rand_end);

      // didn't find an overlap, keep going
      if (overlaps.empty())
        continue;

      // check that the chosen end is <= the end of the overlapping interval
      bool enclosed = true;
      for (auto& j : overlaps) {
        if (rand_start >= j.start) {
          if (rand_end > j.stop) {
            enclosed = false;
          }
        }
      }
      if (!enclosed)
        continue;

      // if we get here, all checks pass. keep the interval.
      inbounds = true;

      starts_out[i] = rand_start;
      ends_out[i] = rand_end;
    }
  }

  return cpp11::writable::data_frame(
      {"chrom"_nm = chroms_out, "start"_nm = starts_out, "end"_nm = ends_out});
}
