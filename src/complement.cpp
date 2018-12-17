// complement.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

//[[Rcpp::export]]
DataFrame complement_impl(ValrGroupedDataFrame gdf, DataFrame genome) {

  genome_map_t chrom_sizes = makeChromSizes(genome) ;

  DataFrame df = gdf.data() ;

  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;
  CharacterVector chroms = df["chrom"] ;

  std::vector<std::string> chroms_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;

  int ngroups = gdf.ngroups() ;
  ListView idx(gdf.indices()) ;

  for (int i = 0; i < ngroups; ++i) {

    IntegerVector indices ;
    indices = idx[i];
    int ni = indices.size() ;

    int start, end ;
    int last_end = 1 ;

    // get chrom from first index
    auto chrom = as<std::string>(chroms[indices[0] - 1]) ;

    for (int j = 0; j < ni; ++j) {

      start = starts[indices[j] - 1] ;
      end = ends[indices[j] - 1] ;

      if (j == 0) {
        if (start == 0) {
          last_end = end ;
          continue ;
        } else {
          chroms_out.push_back(chrom) ;
          starts_out.push_back(0) ;
          ends_out.push_back(start) ;
        }
      } else {
        chroms_out.push_back(chrom) ;
        starts_out.push_back(last_end) ;
        ends_out.push_back(start) ;
      }

      last_end = end;
    }

    auto chrom_size = chrom_sizes[chrom] ;

    if (last_end < chrom_size) {
      chroms_out.push_back(chrom) ;
      starts_out.push_back(last_end) ;
      ends_out.push_back(chrom_size) ;
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
genome <- tibble::tribble(
   ~chrom,  ~size,
   "chr1", 500,
   "chr2", 600,
   "chr3", 800
)

x <- tibble::tribble(
   ~chrom, ~start, ~end,
   "chr1", 100,    300,
   "chr1", 200,    400,
   "chr2", 1,      100,
   "chr2", 200,    400,
   "chr3", 500,    600
) %>% group_by(chrom)

# intervals not covered by x
x <- bed_merge(x) %>% group_by(chrom)

complement_impl(x, genome)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 100
)
bed <- tibble::tribble(
  ~chrom,   ~start,    ~end,
  "chr1",        1,      20
) %>% group_by(chrom)

complement_impl(bed, genome)

*/
