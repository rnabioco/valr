// bed_flank.cpp
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

//[[Rcpp::export]]
DataFrame flank_impl(DataFrame inputTable, DataFrame genome,
                     double both = 0, double left = 0, double right = 0,
                     bool fraction = false, bool strand = false, bool trim = false) {

  // Warnings
  if (both == 0 & left == 0 & right == 0)
    stop("specify one of both, left, right");

  if (both != 0 & (left != 0 || right != 0))
    stop("ambiguous side spec for bed_flank");

  std::vector<std::string> TableNames = inputTable.names();
  int TableLen = TableNames.size();

  bool strandTest = false;

  for (int i = 0; i < TableLen; i++)
    if (TableNames[i] == "strand")
      strandTest = true;

  if (strand == true & strandTest == false)
    stop("expected strand column");


  // Set both
  if (both > 0) left = right = both;


  // Set input and output vectors
  std::vector<std::string> chroms = inputTable["chrom"];
  std::vector<int> startCoords    = inputTable["start"];
  std::vector<int> endCoords      = inputTable["end"];

  int N = startCoords.size();
  std::vector<int> coordSize(N);
  std::vector<int> idxOut;

  std::vector<double> startOut;
  std::vector<double> endOut;


  // Create unordered map for chrom sizes
  genome_map_t chroMap = makeChromSizes(genome);



  for (int i = 0; i < N; i++) {

    int leftstart;
    int leftend;
    int rightstart;
    int rightend;

    // strand
    if (strand == true) {
      std::vector<std::string> strands = inputTable["strand"];

      // strand, fraction
      if (fraction == true) {
        coordSize[i] = endCoords[i] - startCoords[i];

        if (strands[i] == "+") {
          leftstart  = startCoords[i] - coordSize[i] * left;
          leftend    = startCoords[i];
          rightstart = endCoords[i];
          rightend   = endCoords[i] + coordSize[i] * right;

        } else {
          leftstart  = endCoords[i];
          leftend    = endCoords[i] + coordSize[i] * left;
          rightstart = startCoords[i] - coordSize[i] * right;
          rightend   = startCoords[i];
        }

      // strand, no fraction
      } else {

        if (strands[i] == "+") {
          leftstart  = startCoords[i] - left;
          leftend    = startCoords[i];
          rightstart = endCoords[i];
          rightend   = endCoords[i] + right;

        } else {
          leftstart  = endCoords[i];
          leftend    = endCoords[i] + left;
          rightstart = startCoords[i] - right;
          rightend   = startCoords[i];
        }
      }

    // no strand
    } else {

      // no strand, fraction
      if (fraction == true) {
        coordSize[i] = endCoords[i] - startCoords[i];

        leftstart  = startCoords[i] - coordSize[i] * left;
        leftend    = startCoords[i];
        rightstart = endCoords[i];
        rightend   = endCoords[i] + coordSize[i] * right;

      // no strand, no fraction
      } else {
        leftstart  = startCoords[i] - left;
        leftend    = startCoords[i];
        rightstart = endCoords[i];
        rightend   = endCoords[i] + right;
      }
    }


    // Compare new intervals to chrom sizes
    std::string chrom = chroms[i];
    int chrSize = chroMap[chrom];

    if (left > 0 & leftstart > 0 & leftend <= chrSize) {
      startOut.push_back (leftstart);
      endOut.push_back   (leftend);
      idxOut.push_back   (i);

    } else if (trim == true & leftstart > 0 & leftend > chrSize) {
      startOut.push_back (leftstart);
      endOut.push_back   (chrSize);
      idxOut.push_back   (i);

    } else if (trim == true & leftstart <= 0 & leftend <= chrSize) {
      startOut.push_back (1);
      endOut.push_back   (leftend);
      idxOut.push_back   (i);

    } else if (trim == true & leftstart <= 0 & leftend > chrSize) {
      startOut.push_back (1);
      endOut.push_back   (chrSize);
      idxOut.push_back   (i);
    }


    if (right > 0 & rightstart > 0 & rightend <= chrSize) {
      startOut.push_back (rightstart);
      endOut.push_back   (rightend);
      idxOut.push_back   (i);

    } else if (trim == true & rightstart > 0 & rightend > chrSize) {
      startOut.push_back (rightstart);
      endOut.push_back   (chrSize);
      idxOut.push_back   (i);

    } else if (trim == true & rightstart <= 0 & rightend <= chrSize) {
      startOut.push_back (1);
      endOut.push_back   (rightend);
      idxOut.push_back   (i);

    } else if (trim == true & rightstart <= 0 & rightend > chrSize) {
      startOut.push_back (1);
      endOut.push_back   (chrSize);
      idxOut.push_back   (i);
    }
  }


  // Write new DataFrame
  DataFrame outTable = DataFrameSubsetVisitors(inputTable, names(inputTable)).subset(idxOut, "data.frame");

  outTable["start"] = startOut;
  outTable["end"] = endOut;

  return outTable;
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

