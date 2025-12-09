// complement.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <string>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/genome.hpp"

using namespace cpp11::literals;

[[cpp11::register]]
cpp11::writable::data_frame complement_impl(cpp11::data_frame gdf, cpp11::data_frame genome) {
  valr::genome_map_t chrom_sizes = valr::make_genome_map(genome);

  valr::GroupedDataFrame grouped(gdf);
  cpp11::data_frame df = grouped.data();

  cpp11::doubles starts = df["start"];
  cpp11::doubles ends = df["end"];
  cpp11::strings chroms = df["chrom"];

  std::vector<std::string> chroms_out;
  std::vector<double> starts_out;
  std::vector<double> ends_out;

  int ngroups = grouped.ngroups();
  cpp11::list idx = grouped.indices();

  for (int i = 0; i < ngroups; ++i) {
    cpp11::integers indices(idx[i]);
    int ni = indices.size();

    double start, end;
    double last_end = 1;

    // Get chrom from first index (indices are 1-based from R)
    std::string chrom(chroms[indices[0] - 1]);

    for (int j = 0; j < ni; ++j) {
      int row_idx = indices[j] - 1;  // Convert to 0-based
      start = starts[row_idx];
      end = ends[row_idx];

      if (j == 0) {
        if (start == 0) {
          last_end = end;
          continue;
        } else {
          chroms_out.push_back(chrom);
          starts_out.push_back(0);
          ends_out.push_back(start);
        }
      } else {
        chroms_out.push_back(chrom);
        starts_out.push_back(last_end);
        ends_out.push_back(start);
      }

      last_end = end;
    }

    double chrom_size = chrom_sizes[chrom];

    if (last_end < chrom_size) {
      chroms_out.push_back(chrom);
      starts_out.push_back(last_end);
      ends_out.push_back(chrom_size);
    }
  }

  int nrow_out = chroms_out.size();
  cpp11::writable::strings out_chroms(nrow_out);
  cpp11::writable::doubles out_starts(nrow_out);
  cpp11::writable::doubles out_ends(nrow_out);

  for (int i = 0; i < nrow_out; ++i) {
    out_chroms[i] = chroms_out[i];
    out_starts[i] = starts_out[i];
    out_ends[i] = ends_out[i];
  }

  return cpp11::writable::data_frame(
      {"chrom"_nm = out_chroms, "start"_nm = out_starts, "end"_nm = out_ends});
}
