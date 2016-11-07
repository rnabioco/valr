#include "valr.h"

//[[Rcpp::export]]
DataFrame complement_impl(GroupedDataFrame gdf, DataFrame genome) {

  genome_map_t chrom_sizes = makeChromSizes(genome) ;

  DataFrame df = gdf.data() ;

  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;
  CharacterVector chroms = df["chrom"] ;

  std::vector<std::string> chroms_out ;
  std::vector<int> starts_out ;
  std::vector<int> ends_out ;

  int ngroups = gdf.ngroups() ;
  GroupedDataFrame::group_iterator git = gdf.group_begin() ;
  for (int i=0; i<ngroups; ++i, ++git) {

    SlicingIndex indices = *git ;
    int ni = indices.size() ;

    int start, end ;
    int last_end = 1 ;

    // get chrom from first index
    auto chrom = as<std::string>(chroms[indices[0]]) ;

    for (int j=0; j<ni; ++j) {

      start = starts[indices[j]] ;
      end = ends[indices[j]] ;

      if (j == 0) {
        if (start == 1) {
          last_end = end ;
          continue ;
        } else {
          chroms_out.push_back(chrom) ;
          starts_out.push_back(1) ;
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

  return DataFrame::create(Named("chrom") = chroms_out,
                           Named("start") = starts_out,
                           Named("end") = ends_out,
                           Named("stringsAsFactors") = false) ;
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
