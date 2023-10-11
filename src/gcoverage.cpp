// gcoverage.cpp
//
// Copyright (C) 2023 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "valr.h"

typedef std::pair<int, int> posTrack;
typedef std::vector<posTrack> posTracker;

// sort starts and ends  while encoding coverage value
// to sum for overlapping intervals
posTracker collatePositions(const IntegerVector& starts,
                            const IntegerVector& ends) {

  posTracker p;
  auto n = starts.size() ;

  if (n != ends.size()) {
    stop("incompatible start and end vector supplied") ;
  }

  for (int i = 0; i < n; i++) {
    auto s = starts[i];
    auto e = ends[i];
    p.push_back(posTrack(s, 1));
    p.push_back(posTrack(e, -1));
  }

  std::sort(p.begin(), p.end()) ;

  return p;
}

//[[Rcpp::export]]
DataFrame gcoverage_impl(const ValrGroupedDataFrame& gdf,
                         const IntegerVector& max_coords) {

  auto ng = gdf.ngroups() ;
  DataFrame df = gdf.data() ;

  ListView idx(gdf.indices()) ;

  if(max_coords.size() != ng) {
    stop("max_coords must equal the number of groups in data.frame");
  }

  std::vector<int> out_indices, depths, starts, ends;

  for (int i = 0; i < ng; i++) {

    auto max_coord = max_coords[i];

    IntegerVector indices = idx[i];
    indices = indices - 1;

    IntegerVector qstarts = df["start"] ;
    IntegerVector qends   = df["end"] ;
    qstarts = qstarts[indices];
    qends = qends[indices];

    // track single index of a representative interval
    // this will be used to recover grouped column values
    int grp_idx = indices[0];

    auto pos = collatePositions(qstarts, qends);

    int prev_pos = 0;
    int cvg = 0;

    for (auto p:pos) {
      if (p.first > max_coord) {
        warning("Out of bounds interval detected at position: %s \n"
                "  Out of bounds intervals will be ignored",
                p.first);
        break;
      }

      if (p.first == prev_pos) {
        cvg += p.second;
        continue;
      }
      depths.push_back(cvg);
      starts.push_back(prev_pos);
      ends.push_back(p.first);
      out_indices.push_back(grp_idx);

      prev_pos = p.first;
      cvg += p.second;
    }

    if(ends.back() < max_coord) {
      depths.push_back(0);
      starts.push_back(ends.back());
      ends.push_back(max_coord);
      out_indices.push_back(grp_idx);
    }

  }

  // subset original dataframe, note that the grouped column values will be correct
  // but non-grouped column values are no longer matched to each interval and
  // must be dropped on the R side
  DataFrame subset_x = subset_dataframe(df, out_indices) ;

  subset_x["start"] = starts ;
  subset_x["end"] = ends ;
  subset_x[".depth"] = depths ;

  return subset_x ;
}

/*** R
library(dplyr)
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 20, 70, 6, 25, "+",
  "chr1", 50, 100, 1, 25, "-",
  "chr1", 200, 250, 3, 25, "+",
  "chr1", 220, 250, 3, 25, "+",
  "chr2", 80, 130, 5, 25, "-",
  "chr2", 150, 200, 4, 25, "+",
  "chr2", 180, 230, 2, 25, "-",
  "chr2", 190, 230, 2, 25, "-"
) |>  group_by(chrom)

gcoverage_impl(x, max_coords = c(1000, 500)) |> as.data.frame()

*/
