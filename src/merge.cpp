#include "valr.h"

//[[Rcpp::export]]
DataFrame merge_impl(GroupedDataFrame gdf, int max_dist = 0) {

  auto ng = gdf.ngroups() ;

  DataFrame df = gdf.data() ;

  auto nr = df.nrows() ;
  auto nc = df.size() ;

  IntegerVector ids(nr) ;      // store ids
  IntegerVector overlaps(nr) ; // store overlap values

  IntegerVector starts   = df["start"] ;
  IntegerVector ends     = df["end"] ;

  GroupedDataFrame::group_iterator git = gdf.group_begin() ;
  for (int i=0; i<ng; i++, ++git) {

    SlicingIndex indices = *git ;

    intervalVector intervals = makeIntervalVector(df, indices);

    interval_t last_interval = interval_t(0,0,0) ;

    int id, last_id = 0 ; // holds start of the first of intervals to be merged

    for (auto it : intervals) {

      auto idx = it.value ;

      // must be `int` type, not `auto`
      int overlap = intervalOverlap(it, last_interval) ;
      overlaps[idx] = overlap ;

      id = it.start ;

      if (overlap > 0) {
        ids[idx] = last_id ;
      } else if (overlap <= 0 && std::abs(overlap) <= max_dist) {
        ids[idx] = last_id ;
      } else {
        ids[idx] = id ;
        last_id = it.start ;
      }

      last_interval = it ;
    }
  }

  // add two new columns, ids and overlaps
  List out(nc + 2) ;
  CharacterVector onames = df.attr("names") ;
  CharacterVector names(nc + 2) ;

  for (int i=0; i<nc; i++) {
    out[i] = df[i] ;
    names[i] = onames[i] ;
  }

  // add ids
  out[nc] = ids;
  std::string id_name = ".id_merge" ;
  names[nc] = id_name ;

  // add overlaps
  out[nc+1] = overlaps ;
  std::string overlaps_name = ".overlap_merge" ;
  names[nc+1] = overlaps_name;

  out.attr("class") = df.attr("class") ;
  out.attr("row.names") = df.attr("row.names") ;
  out.attr("names") = names ;

  return out ;
}

/*** R
  library(dplyr)
  bed_tbl <- tibble::tribble(
    ~chrom, ~start, ~end, ~value,
    "chr1", 1,      50,   1,
    "chr1", 100,    200,  2,
    "chr1", 150,    250,  3,
    "chr1", 175,    225,  3.5,
    "chr2", 1,      25,   4,
    "chr2", 200,    400,  5,
    "chr2", 400,    500,  6,
    "chr2", 450,    550,  7,
    "chr3", 450,    550,  8,
    "chr3", 500,    600,  9
  ) %>% group_by(chrom)

  merge_impl(bed_tbl) %>% as_data_frame
*/
