// random.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr/random.hpp"

#include <cpp11.hpp>
#include <cpp11/R.hpp>

#include <Rmath.h>

#include <cmath>
#include <numeric>
#include <vector>

using namespace cpp11::literals;

[[cpp11::register]]
cpp11::writable::data_frame random_impl(cpp11::data_frame genome, double length, int n,
                                        int seed = 0) {
  cpp11::strings chroms = genome["chrom"];
  cpp11::doubles sizes = genome["size"];

  int nchrom = chroms.size();

  if (seed == 0) {
    seed = static_cast<int>(std::round(::Rf_runif(0, RAND_MAX)));
  }

  // Seed the generator
  valr::Engine generator(seed);

  // Calculate weights for chrom distribution
  double mass = std::accumulate(sizes.begin(), sizes.end(), 0.0);

  std::vector<double> weights(nchrom);
  for (int i = 0; i < nchrom; ++i) {
    weights[i] = sizes[i] / mass;
  }

  std::vector<int> chromidx(nchrom);
  std::iota(chromidx.begin(), chromidx.end(), 0);

  valr::PiecewiseConstDist chrom_dist(chromidx.begin(), chromidx.end(), weights.begin());

  // Make and store a distribution for each chrom size
  std::vector<valr::UniformIntDist> size_rngs;
  size_rngs.reserve(nchrom);

  for (int i = 0; i < nchrom; ++i) {
    double size = sizes[i];
    // Sub length to avoid off-chrom coordinates
    valr::UniformIntDist size_dist(0, static_cast<int>(size - length));
    size_rngs.push_back(size_dist);
  }

  cpp11::writable::strings rand_chroms(n);
  cpp11::writable::doubles rand_starts(n);

  for (int i = 0; i < n; ++i) {
    int chrom_idx = static_cast<int>(chrom_dist(generator));
    rand_chroms[i] = chroms[chrom_idx];

    valr::UniformIntDist size_dist = size_rngs[chrom_idx];
    double rand_start = size_dist(generator);
    rand_starts[i] = rand_start;
  }

  cpp11::writable::doubles rand_ends(n);
  for (int i = 0; i < n; ++i) {
    rand_ends[i] = rand_starts[i] + length;
  }

  return cpp11::writable::data_frame(
      {"chrom"_nm = rand_chroms, "start"_nm = rand_starts, "end"_nm = rand_ends});
}
