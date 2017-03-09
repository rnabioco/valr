#include "valr.h"

genome_map_t makeChromSizes(DataFrame genome) {

  genome_map_t chrom_sizes ;

  CharacterVector chroms = genome["chrom"] ;
  IntegerVector sizes = genome["size"] ;

  int nchrom = genome.nrows() ;
  for (int i = 0; i < nchrom; ++i) {
    std::string chrom = as<std::string>(chroms[i]) ;
    int size = sizes[i] ;
    chrom_sizes.insert({chrom, size}) ;
  }

  return chrom_sizes ;
}

// the value field of intervals in the returned vector correspond to the index
// of the interval in the original dataframe (i.e., the values of the
// SlicingIndex)
intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) {

  intervalVector intervals ;

  IntegerVector starts = df["start"] ;
  IntegerVector ends   = df["end"] ;

  int size = si.size() ;

  for (int i = 0; i < size; ++i) {
    int j = si[i] ;
    intervals.push_back(interval_t(starts[j], ends[j], j)) ;
  }
  return intervals ;
}

// compare two dataframe rows and report equality
bool compareDataFrameRows(DataFrame df_x, DataFrame df_y, int idx_x, int idx_y) {

  IntegerVector x_idx = IntegerVector::create(idx_x) ;
  IntegerVector y_idx = IntegerVector::create(idx_y) ;

  DataFrame subset_x = DataFrameSubsetVisitors(df_x, names(df_x)).subset(x_idx, "data.frame");
  DataFrame subset_y = DataFrameSubsetVisitors(df_y, names(df_y)).subset(y_idx, "data.frame");

  int ncols = df_x.size() ;
  bool cols_equal = false;
  for (int i = 0; i < ncols; i++) {
    CharacterVector col_x = subset_x[i] ;
    CharacterVector col_y = subset_y[i] ;
    cols_equal = is_true(all(col_x == col_y)) ;
    if (!cols_equal) {
      break;
    }
  }
  return cols_equal ;
}

