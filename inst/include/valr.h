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

// boost icl
#include <boost/icl/interval_set.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>

using namespace boost::icl ;

typedef discrete_interval<int>   icl_interval_t ;
typedef interval_set<int>        icl_interval_set_t ;

extern icl_interval_set_t makeIclIntervalSet(DataFrame df, SlicingIndex si) ;

#endif
