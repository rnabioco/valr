#include "valr.h"

void subtract_group(intervalVector vx, intervalVector vy,
                    std::vector<int>& indices_out,
                    std::vector<int>& starts_out, std::vector<int>& ends_out) {

  intervalTree tree_y(vy) ;
  intervalVector overlaps ;
  IntervalStartSorter<int, int> intervalStartSorter ;

  for (auto it : vx) {

    auto x_start = it.start;
    auto x_stop = it.stop;

    tree_y.findOverlapping(it.start, it.stop, overlaps) ;

    // compute number of overlaps
    int overlap_count = overlaps.size();

    // handle no overlaps and continue
    if (overlap_count == 0) {
      indices_out.push_back(it.value) ;
      starts_out.push_back(it.start) ;
      ends_out.push_back(it.stop) ;
      continue;
    }
    // sort overlaps by start not sure if necessary
    std::sort(overlaps.begin(), overlaps.end(), intervalStartSorter) ;

    // keep track of the number of new intervals to generate
    int new_ivl_count = 0;

    //iterate through overlaps with current x  interval
    // modifying start and stop as necessary
    for (auto oit : overlaps) {

      auto y_start = oit.start;
      auto y_stop = oit.stop;

      if (y_start <= x_start) {
        //advance x_start to end of y
        x_start = std::min(y_stop, x_stop) ;
        continue ;
      } else if (y_start > x_start) {
        // report new interval
        indices_out.push_back(it.value) ;
        starts_out.push_back(x_start) ;
        ends_out.push_back(y_start) ;
        // advance to end of y ivl
        x_start = y_stop;
        new_ivl_count++ ;
        continue ;
      }
    }

    if (x_start < x_stop) {
      indices_out.push_back(it.value) ;
      starts_out.push_back(x_start) ;
      ends_out.push_back(x_stop) ;
    }

    overlaps.clear() ;
  }
}
//[[Rcpp::export]]
DataFrame subtract_impl(GroupedDataFrame gdf_x, GroupedDataFrame gdf_y) {

  std::vector<std::string> chrom_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;

  auto ng_x = gdf_x.ngroups() ;
  auto ng_y = gdf_y.ngroups() ;

  DataFrame df_x = gdf_x.data() ;
  DataFrame df_y = gdf_y.data() ;

  // get labels info for grouping
  DataFrame labels_x(df_x.attr("labels"));
  DataFrame labels_y(df_y.attr("labels"));

  GroupedDataFrame::group_iterator git_x = gdf_x.group_begin() ;
  // indices_to_report
  std::vector<int> indices_out ;
  for (int nx = 0; nx < ng_x; nx++, ++git_x) {

    SlicingIndex indices_x = *git_x ;

    // keep track of if x chrom is present in y
    bool group_seen(false);

    GroupedDataFrame::group_iterator git_y = gdf_y.group_begin() ;
    for (int ny = 0; ny < ng_y; ny++, ++git_y) {

      SlicingIndex indices_y = *git_y ;

      bool same_groups = compareDataFrameRows(labels_x, labels_y, nx, ny);

      if (same_groups) {
        group_seen = true ;

        intervalVector vx = makeIntervalVector(df_x, indices_x) ;
        intervalVector vy = makeIntervalVector(df_y, indices_y) ;

        subtract_group(vx, vy,
                       indices_out,
                       starts_out, ends_out) ;

      }
    }

    // return x intervals if x chromosome not found in y
    if (group_seen) {
      continue;
    } else {
      intervalVector vx = makeIntervalVector(df_x, indices_x) ;
      for (auto it : vx) {
        indices_out.push_back(it.value) ;
        starts_out.push_back(it.start) ;
        ends_out.push_back(it.stop) ;
      }
    }
  }

  // extract out x data, new intervals will be generated as copies of the parent interval
  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(indices_out, "data.frame");

  auto ncol_x = subset_x.size() ;

  CharacterVector names(ncol_x) ;
  CharacterVector names_x = subset_x.attr("names") ;

  // report same number of columns
  List out(ncol_x) ;

  // build new dataframe with colnames and existing data
  for (int i = 0; i < ncol_x; i++) {
    auto name_x = as<std::string>(names_x[i]) ;
    names[i] = name_x ;
    out[i] = subset_x[i] ;
  }
  out.attr("names") = names ;
  out.attr("class") = classes_not_grouped() ;
  auto nrows = subset_x.nrows() ;
  set_rownames(out, nrows) ;

  // assign new starts and ends
  out["start"] = starts_out ;
  out["end"] = ends_out ;

  return out ;

}

/***R
library(dplyr)
library(valr)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 1e7
)

n <- 1e4
x <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)
y <- bed_random(genome, n = n) %>% bed_sort %>% group_by(chrom)

subtract_impl(x, y) %>% as_data_frame()

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    200
) %>% group_by(chrom)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1000,    2000
) %>% group_by(chrom)

subtract_impl(x, y)

*/
