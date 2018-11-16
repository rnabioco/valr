// utils.h
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__utils_H
#define valr__utils_H

#include "valr.h"

DataFrame subset_dataframe(const DataFrame& df,
                           std::vector<int> indices,
                           SEXP frame) ;

DataFrame subset_dataframe(const DataFrame& df,
                           IntegerVector indices,
                           SEXP frame) ;

#endif
