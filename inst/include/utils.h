// utils.h
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__dplyrutils_H
#define valr__dplyrutils_H

#include "valr.h"

DataFrame subset_dataframe(const DataFrame& df,
                           std::vector<int> indices) ;

DataFrame subset_dataframe(const DataFrame& df,
                           IntegerVector indices) ;

inline DataFrame check_is_grouped(const DataFrame& x) {
  bool is_grouped(Rf_inherits(x, "grouped_df")) ;

  if (!is_grouped) {
    Rcpp::stop("error: grouped dataframe required") ;
  }
  return (x) ;
}

template <typename Df>
inline void set_rownames(Df& data, int n) {
  data.attr("row.names") =
    Rcpp::IntegerVector::create(Rcpp::IntegerVector::get_na(), -n);
}

namespace Rcpp {

  typedef Vector<INTSXP, NoProtectStorage> IntegerVectorView;
  typedef Vector<VECSXP, NoProtectStorage> ListView;
  typedef DataFrame_Impl<NoProtectStorage> DataFrameView;

}

#endif
