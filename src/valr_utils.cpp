#include "valr.h"

DataFrame subset_dataframe(const DataFrame& df,
                           std::vector<int> indices,
                           SEXP frame){

  IntegerVector idx_rcpp = wrap(indices) ;
  IntegerVector idx_r = idx_rcpp + 1 ;
  DataFrame out = DataFrameSubsetVisitors(df, frame).subset_all(idx_r);
  return(out) ;
}

DataFrame subset_dataframe(const DataFrame& df,
                           IntegerVector indices,
                           SEXP frame){

  IntegerVector idx_r = indices + 1 ;
  DataFrame out = DataFrameSubsetVisitors(df, frame).subset_all(idx_r);
  return(out) ;
}
