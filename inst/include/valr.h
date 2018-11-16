// valr.h
//
// Copyright (C) 2016 - 2017 Jay Hesselberth and Kent Riemondy
//
// This file is part of valr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#ifndef valr__valr_H
#define valr__valr_H

// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
using namespace Rcpp ;

#include <dplyr/dplyr.h>
#include <dplyr/visitors/subset/DataFrameSubsetVisitors.h>
#include <dplyr/visitors/subset/DataFrameSelect.h>
using namespace dplyr ;

#include "utils.h"
#include "IntervalTree.h"
#include "intervals.h"
#include "group_apply.h"
#include "genome.h"
#include "random.h"
#include "DataFrameBuilder.h"
#include "valr_dplyr_types.h"
#endif
