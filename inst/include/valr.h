#ifndef valr__valr_H 
#define valr__valr_H 

#include <Rcpp.h>
#include <dplyr.h>
#include "IntervalTree.h"
//[[Rcpp::depends(dplyr)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;
using namespace dplyr ;

typedef Interval<int>                interval_t ;
typedef std::vector< interval_t >    intervalVector ;
typedef IntervalTree<int>            intervalTree ;

extern int intervalOverlap(const interval_t& a, const interval_t& b) ;
extern intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) ;

#endif
