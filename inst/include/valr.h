#ifndef valr__valr_H
#define valr__valr_H

#include <Rcpp.h>
#include <dplyr.h>
#include "IntervalTree.h"
#include <functional>
#include <random>
#include <stack>

//[[Rcpp::depends(dplyr)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;
using namespace dplyr ;


typedef std::mt19937                           ENG ;
typedef std::uniform_int_distribution<int>     UDIST ;
typedef std::piecewise_constant_distribution<> PDIST ;

typedef std::unordered_map<std::string, int> genome_map_t ;
extern genome_map_t makeChromSizes(DataFrame genome) ;

typedef Interval<int>                interval_t ;
typedef std::vector< interval_t >    intervalVector ;
typedef IntervalTree<int>            intervalTree ;

extern intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) ;
extern bool compareDataFrameRows(DataFrame x, DataFrame y, int idx_x, int idx_y) ;
// template function to generate matched chromosome
// interval vectors and apply a function to each interval vector

template < typename FN, typename... ARGS >
void PairedGroupApply(const GroupedDataFrame& x,
                      const GroupedDataFrame& y,
                      FN&& fn, ARGS&& ... args) {

  auto data_x = x.data() ;
  auto data_y = y.data() ;

  auto ng_x = x.ngroups() ;
  auto ng_y = y.ngroups() ;

  // get labels info for grouping
  DataFrame labels_x(data_x.attr("labels"));
  DataFrame labels_y(data_y.attr("labels"));

  GroupedDataFrame::group_iterator git_x = x.group_begin() ;
  for (int nx = 0; nx < ng_x; nx++, ++git_x) {

    SlicingIndex gi_x = *git_x ;

    GroupedDataFrame::group_iterator git_y = y.group_begin() ;
    for (int ny = 0; ny < ng_y; ny++, ++git_y) {

      SlicingIndex gi_y = *git_y ;

      // make sure that x and y groups are the same
      bool same_groups = compareDataFrameRows(labels_x, labels_y, nx, ny);

      if (same_groups) {

        intervalVector vx = makeIntervalVector(data_x, gi_x) ;
        intervalVector vy = makeIntervalVector(data_y, gi_y) ;

        std::bind(std::forward<FN>(fn), vx, vy, std::forward<ARGS>(args)...)();
      }
    }
  }
}


#endif
