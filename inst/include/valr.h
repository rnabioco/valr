#ifndef valr__valr_H
#define valr__valr_H

// [[Rcpp::depends(dplyr)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <dplyr.h>
#include <random>
#include <functional>
#include <stack>

#include "genome.h"
#include "intervals.h"
#include "group_apply.h"
#include "IntervalTree.h"

using namespace Rcpp ;
using namespace dplyr ;

typedef std::mt19937                           ENGINE ;
typedef std::uniform_int_distribution<int>     UINT_DIST ;
typedef std::piecewise_constant_distribution<> PCONST_DIST ;

#endif
