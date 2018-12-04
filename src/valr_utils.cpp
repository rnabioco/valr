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

GroupedDataFrame::GroupedDataFrame(DataFrame x):
  data_(check_is_grouped(x)),
  groups(data_.attr("groups"))
{}


DataFrame extract_groups(const DataFrame& x) {
  int n = x.ncol() - 1;
  CharacterVector df_names = x.names() ;

  List res(n);
  CharacterVector new_names(n) ;

  for (int i = 0; i < n; i++) {
    res[i] = x[i] ;
    new_names[i] = df_names[i] ;
  }

  set_rownames(res, x.nrow()) ;
  res.names() = new_names ;
  res.attr("class") = "data.frame" ;
  return res;
}

// Perform intersect between group_data of groupeddataframes and return the
// row indicies shared between both dataframes. Used to identify matching
// groups for two table operations i.e. intersect_impl.
// Based on intersect_data_frame from dplyr
std::vector<int> shared_row_indexes(const GroupedDataFrame& x,
                                    const GroupedDataFrame& y) {

  DataFrame grp_x = extract_groups(x.group_data()) ;
  DataFrame grp_y = extract_groups(y.group_data()) ;

  typedef VisitorSetIndexSet<DataFrameJoinVisitors> Set;

  SymbolVector x_names = grp_x.names();
  DataFrameJoinVisitors visitors(grp_x, grp_y, x_names, x_names, true, true);
  Set set(visitors);

  train_insert(set, grp_x.nrows());

  std::vector<int> indices ;
  int n_y = grp_y.nrows();
  for (int i = 0; i < n_y; i++) {
    Set::iterator it = set.find(-i - 1);
    if (it != set.end()) {
      indices.push_back(*it);
      set.erase(it);
    }
  }
  return indices ;
}


