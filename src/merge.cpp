//
// merge.cpp
//
// algorithm from www.geeksforgeeks.org/merging-intervals/
//

#include <Rcpp.h>
using namespace Rcpp ;

#include <stack>
#include "intervals.h"


void
store_intervals(std::stack<interval_t>& chrom_intervals, std::list<interval_t>& merged_intervals) {
  
  while ( ! chrom_intervals.empty() ) {
    merged_intervals.push_back( chrom_intervals.top() ) ; 
    chrom_intervals.pop() ;
  }
}


std::list <interval_t>
merge_intervals(std::list <interval_t> intervals) {

  interval_t prev_interval = make_interval("", 0, 0) ;
  int interval_count = 0 ;
  int if_case = 0 ;
  
  std::list <interval_t> merged_intervals ;
  std::stack <interval_t> chrom_intervals ;
  
  std::list<interval_t>::iterator it ; 
  for (it = intervals.begin(); it != intervals.end(); ++it) {
 
    interval_t curr_interval = *it ;
   
    Rcerr << "eval interval: " << curr_interval.chrom << ":" << curr_interval.start << "-" << curr_interval.end << std::endl ; 
    // first interval on chrom
    if ( interval_count == 0 ) {
      
      if ( it == intervals.begin() ) {
        Rcerr << "very first interval!" << std::endl ;
        // first interval on first chrom
        chrom_intervals.push(curr_interval) ;
        // ++interval_count ;
        // prev_interval = *it ;
        // continue ;
        // 
      } else {
        // prev_interval is first interval on curr chrom
        Rcerr << "first interval on new chrom" << std::endl ;
        chrom_intervals.push(prev_interval) ;
      }
      
      if_case = 1 ;   
    }
   
    // switched chromosomes
    else if ( curr_interval.chrom != prev_interval.chrom ) {
      
      Rcerr << "switched chroms" << std::endl ;
      // store current set of intervals
      store_intervals(chrom_intervals, merged_intervals) ;
      chrom_intervals.push(curr_interval) ;
      interval_count = 0 ;
      // prev_interval = *it ;
      // continue ;
      // 
      if_case = 2 ;   
    } 

    else if (interval_overlap(curr_interval, chrom_intervals.top()) > 0) {
      // update the stack interval with new end
      chrom_intervals.top().end = curr_interval.end ;
      
      if_case = 3 ;   
    }
    
    else {
      chrom_intervals.push(curr_interval) ;
      if_case = 4 ;   
    }
    
    Rcerr << "if case: " << if_case << std::endl ;
    
    ++interval_count ;
    prev_interval = *it ;
    
  }
 
  // store last set of intervals 
  store_intervals(chrom_intervals, merged_intervals) ;

  return merged_intervals ;
}


// [[Rcpp::export]]
Rcpp::DataFrame
merge_impl(DataFrame df) {

  // should be replaced with efficient iterator  
  std::list <interval_t> intervals = create_intervals(df) ;

  Rcpp::CharacterVector chroms_v ;
  Rcpp::NumericVector starts_v ;    
  Rcpp::NumericVector ends_v ;    

  std::list<interval_t> merged_intervals = merge_intervals(intervals) ;
 
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
