// DataFrameBuilder.h
//
// Copyright (C) 2016 - 2022 Jay Hesselberth and Kent Riemondy
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
  List data ;
  DataFrameBuilder() {} ;

  // output List with:  List out = DataFrameBuilder
  inline operator List() const {
    return data ;
  }

  inline size_t size() const {
    return data.size() ;
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
      this->names.push_back(name) ;
      this->data.push_back(df[i]) ;
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
      this->names.push_back(name) ;
      this->data.push_back(df[i]) ;
    }
  }

  // apply common attributes to output dataframe
  // by default groups are stripped and tbl_df is returned
  inline List format_df(int nrow) {
    auto names = this->names ;

    set_rownames(this->data, nrow) ;
    this->data.attr("names") = names ;

    if (Rf_inherits(this->data, "grouped_df")) {
      ValrGroupedDataFrame::strip_groups(this->data) ;
    }

    this->data.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame") ;
    return this->data ;
  }
};

#endif
