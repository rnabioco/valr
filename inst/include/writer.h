// writer.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__writer_H
#define valr__writer_H

#include "valr.h"

class DataFrameBuilder {
public:
  std::vector<std::string> names ;
  std::vector<SEXP> data ; // set to SEXP so that it can handle any R type vector
  DataFrameBuilder(){} ;

  // non-SEXP objects should be pass with wrap(obj)
  inline void add_vec(std::string name, SEXP x){
    names.push_back(name) ;
    data.push_back(x) ;
  }

  inline operator List() const {
    List out = wrap(data) ;
    return out ;
  }

  inline size_t size() const {
    return data.size() ;
  }

  inline void add_df(const DataFrame& df,
                     std::string suffix,
                     bool drop_chrom = true){

    auto ncol = df.size() ;

    //  CharacterVector names
    CharacterVector names = df.attr("names") ;

    // names, data
    for (int i = 0; i < ncol; i++) {
      auto name = as<std::string>(names[i]) ;
      if (name != "chrom") {
        name += suffix ;
      } else if (drop_chrom){
        continue ;
      }

      this->add_vec(name, df[i]) ;
    }
  }
};

inline List structure_df(const DataFrameBuilder& out,
                  int nrow){
  List res = out ;
  auto names = out.names ;
  set_rownames(res, nrow ) ;
  res.attr("names") = names ;
  res.attr("class") = classes_not_grouped() ;
  return res ;
}



#endif
