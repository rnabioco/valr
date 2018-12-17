#include "valr.h"

// [[Rcpp::export]]
// Based on Kevin Ushey's implementation here
// http://kevinushey.github.io/blog/2015/01/24/understanding-data-frame-subsetting/
DataFrame rowwise_subset_df(const DataFrame& x,
                            IntegerVector row_indices) {

  // Get some useful variables
  // (lengths of index vectors)
  int column_indices_n = x.ncol();
  int row_indices_n = row_indices.size();

  // Translate from R to C indices.
  // This could be optimized...
  row_indices = row_indices - 1;

  // Allocate the output 'data.frame'.
  List output = no_init(column_indices_n);

  // Set the names, based on the names of 'x'.
  // We'll assume it has names.
  CharacterVector x_names =
    as<CharacterVector>(Rf_getAttrib(x, R_NamesSymbol));

  // Use Rcpp's sugar subsetting to subset names.
  output.attr("names") = x_names;

  // Loop and fill!
  for (int j = 0; j < column_indices_n; ++j)
  {
    // Get the j'th element of 'x'. We don't need to
    // PROTECT it since it's implicitly protected as a
    // child of 'x'.
    SEXP element = VECTOR_ELT(x, j);

    // Get the 'rows' for that vector, and fill.
    SEXP vec = PROTECT(
      Rf_allocVector(TYPEOF(element), row_indices_n)
    );

    for (int i = 0; i < row_indices_n; ++i)
    {
      // Copying vectors is a pain in the butt, because
      // we need to know the actual type underneath the
      // SEXP. Raw and Complex vectors are not handled.
      // Also, if the row_indices contain NA, return type dependent NAs
      switch (TYPEOF(vec))
      {
      case REALSXP:
        if(row_indices[i] == NA_INTEGER) {
          REAL(vec)[i] = Rcpp::Vector<REALSXP>::get_na() ;
        } else {
          REAL(vec)[i] =
            REAL(element)[row_indices[i]];
        }
        break;
      case INTSXP:
      case LGLSXP:
        if(row_indices[i] == NA_INTEGER) {
          INTEGER(vec)[i] = Rcpp::Vector<INTSXP>::get_na() ;
        } else {
          INTEGER(vec)[i] =
            INTEGER(element)[row_indices[i]];
        }
        break;
      case STRSXP:
        if(row_indices[i] == NA_INTEGER) {
          SET_STRING_ELT(vec, i,
                         Rcpp::Vector<STRSXP>::get_na()) ;
          break;
        }
        SET_STRING_ELT(vec, i,
                       STRING_ELT(element, row_indices[i]));
        break;
      case VECSXP:
        if(row_indices[i] == NA_INTEGER) {
          SET_VECTOR_ELT(vec, i,
                         Rcpp::Vector<VECSXP>::get_na()) ;
        } else {
          SET_VECTOR_ELT(vec, i,
                         VECTOR_ELT(element, row_indices[i]));
        }
        break;
      default: {
          stop("Incompatible column type detected");
        }
      }
    }

    // Don't need to protect 'vec' anymore
    UNPROTECT(1);

    // And make sure the output list now
    // refers to that vector!
    SET_VECTOR_ELT(output, j, vec);
  }

  // Finally, copy the attributes of `x` to `output`...
  Rf_copyMostAttrib(x, output);

  // ... but set the row names manually. Note that this
  // is the secret method for creating 'compact' row
  // names, whereby internally R stores the 'empty' row
  // names object as 'c(NA_integer_, -<nrow>)'.
  Rf_setAttrib(output, R_RowNamesSymbol,
               IntegerVector::create(NA_INTEGER, -row_indices_n));

  return output;

}

DataFrame subset_dataframe(const DataFrame& df,
                           std::vector<int> indices,
                           SEXP frame){

  IntegerVector idx_rcpp = wrap(indices) ;
  IntegerVector idx_r = idx_rcpp + 1 ;
  DataFrame out = rowwise_subset_df(df, idx_r);
  return(out) ;
}

DataFrame subset_dataframe(const DataFrame& df,
                           IntegerVector indices,
                           SEXP frame){

  IntegerVector idx_r = indices + 1 ;
  DataFrame out = rowwise_subset_df(df, idx_r);
  return(out) ;
}

// GroupedDataFrame::GroupedDataFrame(DataFrame x):
//   data_(check_is_grouped(x)),
//   groups(data_.attr("groups"))
// {}


ValrGroupedDataFrame::ValrGroupedDataFrame(DataFrame x):
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
// std::vector<int> shared_row_indexes(const ValrGroupedDataFrame& x,
//                                     const ValrGroupedDataFrame& y) {
//
//   DataFrame grp_x = extract_groups(x.group_data()) ;
//   DataFrame grp_y = extract_groups(y.group_data()) ;
//
//   typedef VisitorSetIndexSet<DataFrameJoinVisitors> Set;
//
//   SymbolVector x_names = grp_x.names();
//   DataFrameJoinVisitors visitors(grp_x, grp_y, x_names, x_names, true, true);
//   Set set(visitors);
//
//   train_insert(set, grp_x.nrows());
//
//   std::vector<int> indices ;
//   int n_y = grp_y.nrows();
//   for (int i = 0; i < n_y; i++) {
//     Set::iterator it = set.find(-i - 1);
//     if (it != set.end()) {
//       indices.push_back(*it);
//       set.erase(it);
//     }
//   }
//   return indices ;
// }


