// valr.h
//
// Copyright (C) 2016 - 2025 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__valr_H
#define valr__valr_H

// cpp11 headers (for new cpp11 code)
#include <cpp11.hpp>
#include <cpp11/R.hpp>
// Note: Do NOT use 'using namespace cpp11' to avoid conflicts with Rcpp

// Rcpp headers (temporary, for existing code during migration)
#include <Rcpp.h>
using namespace Rcpp;

// Standard library
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>

// valr headers (old, will be replaced)
#include "utils.h"
#include "grouped_dataframe.h"
#include "IntervalTree.h"
#include "intervals.h"
#include "IntervalTree_ext.h"
#include "group_apply.h"
#include "genome.h"
#include "random.h"
#include "DataFrameBuilder.h"

#endif
