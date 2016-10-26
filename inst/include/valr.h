#ifndef valr__valr_H 
#define valr__valr_H 

#include <Rcpp.h>
#include <dplyr.h>
#include "IntervalTree.h"
#include <functional>
#include <random>

// [[Rcpp::depends(dplyr)]]

using namespace Rcpp ;
using namespace dplyr ;

#include "boost/random.hpp"

typedef boost::random::mt19937                           ENG ;
typedef boost::random::uniform_int_distribution<int>     UDIST ;
typedef boost::random::piecewise_constant_distribution<> PDIST ;

typedef std::map<std::string, int> genome_map_t ;
extern genome_map_t makeChromSizes(DataFrame genome) ;

typedef Interval<int, int>           interval_t ;
typedef std::vector< interval_t >    intervalVector ;
typedef IntervalTree<int>            intervalTree ;

extern intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) ;
extern bool compareDataFrameRows(DataFrame x, DataFrame y, int idx_x, int idx_y ) ;

#endif
