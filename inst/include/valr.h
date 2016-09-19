#ifndef valr__valr_H 
#define valr__valr_H 

#include <Rcpp.h>
#include <dplyr.h>
#include "IntervalTree.h"
#include <functional>
#include <random>
#include <boost/icl/interval_set.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>

//[[Rcpp::depends(dplyr)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;
using namespace dplyr ;
using namespace boost::icl ;

typedef std::mt19937                           ENG ;
typedef std::uniform_int_distribution<int>     UDIST ;
typedef std::piecewise_constant_distribution<> PDIST ;

typedef std::unordered_map<std::string, int> genome_map_t ;
extern genome_map_t makeChromSizes(DataFrame genome) ;

typedef Interval<int>                interval_t ;
typedef std::vector< interval_t >    intervalVector ;
typedef IntervalTree<int>            intervalTree ;

extern intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) ;

// template function to generate matched chromosome
// interval vectors and apply a function to each interval vector

template < typename FN, typename... ARGS >
void PairedGroupApply( const GroupedDataFrame& x, 
                const GroupedDataFrame& y,
                FN&& fn, ARGS&&... args ){
  
  auto data_x = x.data() ;
  auto data_y = y.data() ;

  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;

  CharacterVector chrom_x = data_x["chrom"];
  CharacterVector chrom_y = data_y["chrom"];
  
    
  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for( int nx=0; nx<ng_x; nx++, ++git_x){
  
    SlicingIndex gi_x = *git_x ;
    int group_x = gi_x[0];
  
    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for( int ny=0; ny<ng_y; ny++, ++git_y) {
    
      SlicingIndex gi_y = *git_y ;
      int group_y = gi_y[0];
    
      if( chrom_x[group_x] == chrom_y[group_y] ) {
      
        intervalVector vx = makeIntervalVector(data_x, gi_x) ;
        intervalVector vy = makeIntervalVector(data_y, gi_y) ;
      
        std::bind( std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)... )();
      }
    } 
  }
}

// boost icl
typedef discrete_interval<int>   icl_interval_t ;
typedef interval_set<int>        icl_interval_set_t ;

extern icl_interval_set_t makeIclIntervalSet(DataFrame df, SlicingIndex si) ;

#endif
