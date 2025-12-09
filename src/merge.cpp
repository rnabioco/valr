// merge.cpp
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

cpp11::writable::data_frame collapseMergedIntervals(cpp11::data_frame gdf, int max_dist) {
  valr::GroupedDataFrame grouped(gdf);
  cpp11::data_frame df = grouped.data();

  int ng = grouped.ngroups();
  cpp11::list idx = grouped.indices();

  // approach from http://www.geeksforgeeks.org/merging-intervals/
  std::vector<valr::ivl_t> s;

  for (int i = 0; i < ng; i++) {
    cpp11::integers indices(idx[i]);

    valr::ivl_vector_t intervals = valr::makeIntervalVector(df, indices);

    // set first interval
    s.push_back(intervals[0]);
    intervals.erase(intervals.begin());

    for (auto& it : intervals) {
      auto top = s.back();
      if (top.stop + max_dist < it.start) {
        // no overlap push to stack
        s.push_back(it);
      } else if (top.stop < it.stop) {
        // overlaps and need to update stack top position
        top.stop = it.stop;
        s.pop_back();
        s.push_back(top);
      }
    }
  }

  std::vector<int> indices_x;
  std::vector<double> group_starts;
  std::vector<double> group_ends;

  // iterate through vector of merged intervals and write to dataframe
  for (auto& it : s) {
    indices_x.push_back(it.value);
    group_starts.push_back(it.start);
    group_ends.push_back(it.stop);
  }

  cpp11::writable::data_frame subset_x = valr::subset_dataframe(df, indices_x);

  // Get column names to find start/end positions
  cpp11::strings col_names(subset_x.attr("names"));
  int ncol = subset_x.size();
  int nrow_out = group_starts.size();

  // Build output vectors
  cpp11::writable::doubles new_starts(nrow_out);
  cpp11::writable::doubles new_ends(nrow_out);

  for (int i = 0; i < nrow_out; i++) {
    new_starts[i] = group_starts[i];
    new_ends[i] = group_ends[i];
  }

  // Build output list, replacing start/end columns
  cpp11::writable::list out_list(ncol);
  cpp11::writable::strings out_names(ncol);

  for (int j = 0; j < ncol; j++) {
    std::string name(col_names[j]);
    out_names[j] = name;

    if (name == "start") {
      out_list[j] = new_starts;
    } else if (name == "end") {
      out_list[j] = new_ends;
    } else {
      out_list[j] = subset_x[j];
    }
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}

cpp11::writable::data_frame clusterMergedIntervals(cpp11::data_frame gdf, int max_dist) {
  valr::GroupedDataFrame grouped(gdf);
  cpp11::data_frame df = grouped.data();

  int ng = grouped.ngroups();
  int nr = grouped.nrows();

  cpp11::writable::integers ids(nr);
  int cluster_id = 0;

  cpp11::list idx = grouped.indices();

  for (int i = 0; i < ng; i++) {
    cpp11::integers indices(idx[i]);
    int ni = indices.size();

    valr::ivl_vector_t intervals = valr::makeIntervalVector(df, indices);

    // store an interval to ensure first interval maintained
    valr::ivl_t last_interval(0, 0, 0);
    valr::ivl_t top = last_interval;

    for (int j = 0; j < ni; j++) {
      valr::ivl_t it = intervals[j];
      // get index to store cluster ids in vector
      int row_idx = it.value;
      last_interval = it;  // set interval to compare

      if (j == 0 || top.stop + max_dist < it.start) {
        // no overlap, update interval and cluster id
        top = it;
        cluster_id++;
        ids[row_idx] = cluster_id;
      } else {
        // overlaps or contained in stack top ivl
        if (it.stop > top.stop) {
          top.stop = it.stop;  // update end position
        }
        ids[row_idx] = cluster_id;
      }
    }
  }

  // Build output using DataFrameBuilder
  valr::DataFrameBuilder out;

  // Add all columns from df (don't drop chrom)
  out.add_df(df, "", false);

  // Add ids column
  out.names.push_back(".id_merge");
  out.data.push_back(ids);

  return out.build(nr);
}

[[cpp11::register]]
cpp11::writable::data_frame merge_impl(cpp11::data_frame gdf, int max_dist = 0,
                                       bool collapse = true) {
  if (!collapse) {
    // return a cluster id per input interval
    return clusterMergedIntervals(gdf, max_dist);
  } else {
    // return only merged intervals
    return collapseMergedIntervals(gdf, max_dist);
  }
}
