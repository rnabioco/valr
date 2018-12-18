// random.h
//
// Copyright (C) 2016 - 2018 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__random_H
#define valr__random_H

#include <random>

typedef std::mt19937                           ENGINE ;
typedef std::uniform_int_distribution<int>     UINT_DIST ;
typedef std::piecewise_constant_distribution<> PCONST_DIST ;

#endif
