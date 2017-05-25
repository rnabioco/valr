// flank.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

//[[Rcpp::export]]
DataFrame flank_impl(DataFrame df, DataFrame genome,
                     double both = 0, double left = 0, double right = 0,
                     bool fraction = false, bool strand = false, bool trim = false) {

  // Warnings
  if (both == 0 && left == 0 && right == 0)
    stop("specify one of both, left, right");

  if (both != 0 && (left != 0 || right != 0))
    stop("ambiguous side spec for bed_flank");

  std::vector<std::string> df_names = df.names();
  bool stranded = false;

  for (int i = 0; i < df_names.size(); i++)
    if (df_names[i] == "strand")
      stranded = true;

  if (strand == true && stranded == false)
    stop("expected strand column");

  // Set both
  if (both > 0) left = right = both;

  // Set input and output vectors
  std::vector<std::string> chroms = df["chrom"];
  std::vector<int> starts = df["start"];
  std::vector<int> ends = df["end"];

  std::vector<int> df_idx;

  std::vector<int> starts_out;
  std::vector<int> ends_out;

  // Create unordered map for chrom sizes
  genome_map_t chrom_sizes = makeChromSizes(genome);

  for (int i = 0; i < starts.size(); i++) {

    int leftstart, leftend, rightstart, rightend ;

    int start = starts[i] ;
    int end = ends[i] ;
    int size = end - start;

    // strand
    if (strand == true) {
      std::vector<std::string> strands = df["strand"];

      // strand, fraction
      if (fraction == true) {

        if (strands[i] == "+") {
          leftstart  = start - size * left;
          leftend    = start;
          rightstart = end;
          rightend   = end + size * right;

        } else {
          leftstart  = end;
          leftend    = end + size * left ;
          rightstart = start - size * right ;
          rightend   = start ;
        }

        // strand, no fraction
      } else {

        if (strands[i] == "+") {
          leftstart  = start - left;
          leftend    = start;
          rightstart = end;
          rightend   = end + right;

        } else {
          leftstart  = end;
          leftend    = end + left;
          rightstart = start - right;
          rightend   = start;
        }
      }

      // no strand
    } else {

      // no strand, fraction
      if (fraction == true) {

        leftstart  = start - size * left;
        leftend    = start;
        rightstart = end;
        rightend   = end + size * right;

        // no strand, no fraction
      } else {
        leftstart  = start - left;
        leftend    = start;
        rightstart = end;
        rightend   = end + right;
      }
    }

    // Compare new intervals to chrom sizes
    std::string chrom = chroms[i];
    int chrom_size = chrom_sizes[chrom];

    if (left > 0 && leftstart > 0 && leftend <= chrom_size) {
      starts_out.push_back(leftstart);
      ends_out.push_back(leftend);
      df_idx.push_back(i);

    } else if (trim == true && leftstart > 0 && leftend > chrom_size) {
      starts_out.push_back(leftstart);
      ends_out.push_back(chrom_size);
      df_idx.push_back(i);

    } else if (trim == true && leftstart <= 0 && leftend <= chrom_size) {
      starts_out.push_back(1);
      ends_out.push_back(leftend);
      df_idx.push_back(i);

    } else if (trim == true && leftstart <= 0 && leftend > chrom_size) {
      starts_out.push_back(1);
      ends_out.push_back(chrom_size);
      df_idx.push_back(i);
    }

    if (right > 0 && rightstart > 0 && rightend <= chrom_size) {
      starts_out.push_back(rightstart);
      ends_out.push_back(rightend);
      df_idx.push_back(i);

    } else if (trim == true && rightstart > 0 && rightend > chrom_size) {
      starts_out.push_back(rightstart);
      ends_out.push_back(chrom_size);
      df_idx.push_back(i);

    } else if (trim == true && rightstart <= 0 && rightend <= chrom_size) {
      starts_out.push_back(1);
      ends_out.push_back(rightend);
      df_idx.push_back(i);

    } else if (trim == true && rightstart <= 0 && rightend > chrom_size) {
      starts_out.push_back(1);
      ends_out.push_back(chrom_size);
      df_idx.push_back(i);
    }
  }

  // Write new DataFrame
  DataFrame out = DataFrameSubsetVisitors(df, names(df)).subset(df_idx, "data.frame");

  out["start"] = starts_out;
  out["end"] = ends_out;

  return out;
}



/*** R
setwd('~/Downloads/')
library(valr)
library(dplyr)

genome <- read_genome('~/Documents/GitHub/valr/inst/extdata/genome.txt.gz')
x <- bed_random(n = 10000000, genome)

x <- read_bed('1k_hg19_gene_names.bed', n_fields = 6)

system.time(bed_flank(x, genome, right = 1000, left = 5000, strand = F, trim = T, fraction = T))

y <- bed_flank(x, genome, right = 300, left = 200, strand = T, trim = T, fraction = F)

*/

