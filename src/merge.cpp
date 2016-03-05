//
// merge.cpp
//
// algorithm from www.geeksforgeeks.org/merging-intervals/
//

#include <Rcpp.h>
#include <stack>
#include "Rbedtools.h"

using namespace Rcpp ;

void
store_intervals(std::stack<interval_t>& chrom_intervals, std::list<interval_t>& merged_intervals) {
  
  while ( ! chrom_intervals.empty() ) {
    merged_intervals.push_back( chrom_intervals.top() ) ; 
    chrom_intervals.pop() ;
  }
}


std::list <interval_t>
merge_intervals(std::list <interval_t> intervals, int max_dist) {

  // make an empty previous to start with
  interval_t previous = make_interval("", 0, 0) ;
  
  std::list <interval_t> merged_intervals ;
  std::stack <interval_t> chrom_intervals ;
  
  int overlap = 0 ;
  
  std::list<interval_t>::iterator it ; 
  for (it = intervals.begin(); it != intervals.end(); ++it) {
 
    interval_t current = *it ;
    
    if ( ! chrom_intervals.empty() ) {
      overlap = interval_overlap(current, chrom_intervals.top() ) ;
    }
  
    // switch chroms, can "switch" onto first chrom from null first previous
    if ( current.chrom != previous.chrom ) {
      
      // store current set of intervals
      store_intervals(chrom_intervals, merged_intervals) ;
      
      // store curr, which is the first on the new chrom
      chrom_intervals.push(current) ;
    } 

    else if ( ! isnan(overlap) && (overlap > 0 || std::abs(overlap) < max_dist )) {
      // update the stack interval with new end
      chrom_intervals.top().end = current.end ;
    }
    
    else {
      chrom_intervals.push(current) ;
    }
    
    previous = *it ;
    
  }
 
  // store last set of intervals 
  store_intervals(chrom_intervals, merged_intervals) ;

  return merged_intervals ;
}


// [[Rcpp::export]]
Rcpp::DataFrame
merge_impl(DataFrame df, int max_dist) {

  // should be replaced with efficient iterator  
  std::list <interval_t> intervals = create_intervals(df) ;

  Rcpp::CharacterVector chroms_v ;
  Rcpp::NumericVector starts_v ;    
  Rcpp::NumericVector ends_v ;    

  std::list<interval_t> merged_intervals = merge_intervals(intervals, max_dist) ;
 
  std::list<interval_t>::const_iterator it; 
  for (it = merged_intervals.begin(); it != merged_intervals.end(); ++it) {
   
    interval_t interval = *it;
  
    chroms_v.push_back(interval.chrom) ; 
    starts_v.push_back(interval.start) ;
    ends_v.push_back(interval.end) ;
   
 }
 
 return DataFrame::create( Named("chrom") = chroms_v,
                           Named("start") = starts_v,
                           Named("end") = ends_v) ;
 
}
