// partition.cpp
//
// Copyright (C) 2018 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include <cpp11.hpp>

#include <algorithm>
#include <vector>

#include "valr/dataframe.hpp"
#include "valr/intervals.hpp"

using namespace cpp11::literals;

class IntervalCache {
  // store vector of intervals to partition
  // and max stop of the stored intervals
 public:
  valr::ivl_vector_t ivls;
  int max_stop = 0;
  void clear() {
    ivls.clear();
    max_stop = 0;
  }
  void push_back(valr::ivl_t ivl) { ivls.push_back(ivl); }
  size_t size() { return ivls.size(); }
};

void partitionIntervals(const IntervalCache& ivl_cache, valr::ivl_vector_t& ivl_result) {
  // convert interval to points
  std::vector<int> ivl_points;
  for (const auto& it : ivl_cache.ivls) {
    ivl_points.push_back(it.start);
    ivl_points.push_back(it.stop);
  }

  // get unique positions
  std::sort(ivl_points.begin(), ivl_points.end());
  ivl_points.erase(std::unique(ivl_points.begin(), ivl_points.end()), ivl_points.end());

  // generate partitions as intervals between
  size_t n_pts = ivl_points.size();
  for (size_t i = 0; i < n_pts - 1; i++) {
    valr::ivl_t new_ivl(ivl_points[i], ivl_points[i + 1], ivl_cache.ivls[0].value);
    ivl_result.push_back(new_ivl);
  }
}

[[cpp11::register]]
cpp11::writable::data_frame partition_impl(cpp11::data_frame gdf, int max_dist = -1) {
  valr::GroupedDataFrame grouped(gdf);
  cpp11::data_frame df = grouped.data();

  int ng = grouped.ngroups();
  cpp11::list idx = grouped.indices();

  std::vector<valr::ivl_t> out_ivls;
  IntervalCache ivl_cache;

  for (int i = 0; i < ng; i++) {
    cpp11::integers indices(idx[i]);
    valr::ivl_vector_t intervals = valr::makeIntervalVector(df, indices);

    // set first interval
    ivl_cache.push_back(intervals[0]);
    ivl_cache.max_stop = intervals[0].stop;
    intervals.erase(intervals.begin());

    for (const auto& it : intervals) {
      auto max_stop = ivl_cache.max_stop;
      if (max_stop + max_dist < it.start) {
        // doesn't overlap
        // clear out cache
        if (ivl_cache.ivls.size() > 1) {
          partitionIntervals(ivl_cache, out_ivls);
          ivl_cache.clear();
        } else {
          out_ivls.push_back(ivl_cache.ivls[0]);
          ivl_cache.clear();
        }

        ivl_cache.push_back(it);
        ivl_cache.max_stop = it.stop;
      }

      else if (max_stop + max_dist < it.stop) {
        // overlaps
        // update max_stop position and store ivl
        ivl_cache.max_stop = it.stop;
        ivl_cache.push_back(it);
      }

      else {
        // contained interval
        // store interval
        ivl_cache.push_back(it);
      }
    }

    // write out final stored interval set
    if (ivl_cache.size() > 1) {
      partitionIntervals(ivl_cache, out_ivls);
      ivl_cache.clear();
    } else {
      out_ivls.push_back(ivl_cache.ivls[0]);
      ivl_cache.clear();
    }
  }

  std::vector<int> indices_x;
  std::vector<double> group_starts;
  std::vector<double> group_ends;

  // iterate through vector of partitioned intervals and write to dataframe
  for (const auto& it : out_ivls) {
    indices_x.push_back(it.value);
    group_starts.push_back(it.start);
    group_ends.push_back(it.stop);
  }

  // subset original dataframe, note that the grouped column values will be correct
  // but non-grouped column values are no longer matched to each interval
  // and are dropped on the R side
  cpp11::writable::data_frame out = valr::subset_dataframe(df, indices_x);

  // Get column names to find start/end positions
  cpp11::strings col_names(out.attr("names"));
  int ncol = out.size();
  int nrow_out = group_starts.size();

  // Build new start/end vectors
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
      out_list[j] = out[j];
    }
  }

  // Set attributes for data frame
  out_list.attr("names") = out_names;
  out_list.attr("class") = "data.frame";
  out_list.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -nrow_out});

  return cpp11::writable::data_frame(out_list);
}
