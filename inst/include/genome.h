// genome.h
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__genome_H
#define valr__genome_H

#include "valr.h"

typedef std::unordered_map<std::string, int> genome_map_t ;

inline genome_map_t makeChromSizes(DataFrame genome,
                                   std::string col_chrom = "chrom",
                                   std::string col_size = "size") {

  genome_map_t chrom_sizes ;

  CharacterVector refs = genome[col_chrom] ;
  IntegerVector sizes = genome[col_size] ;

  if (unique(refs).length() != refs.length())
    stop("duplicate reference names in genome file.") ;

  int nchrom = genome.nrows() ;
  for (int i = 0; i < nchrom; ++i) {
    std::string ref = as<std::string>(refs[i]) ;
    int size = sizes[i] ;
    chrom_sizes.insert({ref, size}) ;
  }

  return chrom_sizes ;
}

#endif
