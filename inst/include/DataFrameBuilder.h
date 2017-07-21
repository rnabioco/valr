// DataFrameBuilder.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__DataFrameBuilder_H
#define valr__DataFrameBuilder_H

#include "valr.h"

class DataFrameBuilder {
public:
  std::vector<std::string> names ;
  std::vector<SEXP> data ; // set to SEXP to handle any R type vector
  DataFrameBuilder() {} ;

  // output List with:  List out = DataFrameBuilder
  inline operator List() const {
    List out = wrap(data) ;
    return out ;
  }

  inline size_t size() const {
    return data.size() ;
  }

  // add vector to DataFrameBuilder
  // non-SEXP objects should be passed with Rcpp::wrap(obj)
  inline void add_vec(std::string name, SEXP x) {
    names.push_back(name) ;
    data.push_back(x) ;
  }

  // add dataframe to DataFrameBuilder with suffix
  inline void add_df(const DataFrame& df,
                     std::string suffix = "",
                     bool drop_chrom = true) {

    auto ncol = df.size() ;
    CharacterVector names = df.attr("names") ;

    for (int i = 0; i < ncol; i++) {
      auto name = as<std::string>(names[i]) ;
      if (name != "chrom") {
        name += suffix ;
      } else if (drop_chrom) {
        continue ;
      }

      this->add_vec(name, df[i]) ;
    }
  }

  // add dataframe to DataFrameBuilder without suffix
  // note overloading add_df definition
  inline void add_df(const DataFrame& df,
                     bool drop_chrom = true) {

    auto ncol = df.size() ;

    CharacterVector names = df.attr("names") ;

    for (int i = 0; i < ncol; i++) {
      auto name = as<std::string>(names[i]) ;
      if (name == "chrom" && drop_chrom) {
        continue ;
      }
      this->add_vec(name, df[i]) ;
    }
  }

  // apply common  attributes to output dataframe
  inline List format_df(int nrow) {
    List res = *this ;
    auto names = this->names ;
    set_rownames(res, nrow) ;
    res.attr("names") = names ;
    res.attr("class") = classes_not_grouped() ;
    return res ;
  }
};

#endif
