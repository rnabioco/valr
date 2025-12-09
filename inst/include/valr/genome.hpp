// genome.hpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef VALR_GENOME_HPP
#define VALR_GENOME_HPP

#include <cpp11.hpp>

#include <string>
#include <unordered_map>
#include <unordered_set>

namespace valr {

// Map from chromosome name to size (double for large genomes)
using genome_map_t = std::unordered_map<std::string, double>;

// Build a genome map from a data frame with chrom and size columns
inline genome_map_t make_genome_map(const cpp11::data_frame& genome,
                                    const std::string& col_chrom = "chrom",
                                    const std::string& col_size = "size") {
  genome_map_t chrom_sizes;

  cpp11::strings refs = genome[col_chrom.c_str()];
  cpp11::doubles sizes = genome[col_size.c_str()];

  int nchrom = refs.size();

  // Check for duplicates
  std::unordered_set<std::string> seen;
  for (int i = 0; i < nchrom; ++i) {
    std::string ref(refs[i]);
    if (seen.count(ref)) {
      cpp11::stop("duplicate reference names in genome file.");
    }
    seen.insert(ref);
  }

  // Build the map
  for (int i = 0; i < nchrom; ++i) {
    std::string ref(refs[i]);
    double size = sizes[i];
    chrom_sizes.emplace(ref, size);
  }

  return chrom_sizes;
}

}  // namespace valr

// Legacy typedef for compatibility during migration
using genome_map_t = valr::genome_map_t;

#endif  // VALR_GENOME_HPP
