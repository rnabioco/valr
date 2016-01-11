//
// intersect.cpp
//
// implements the "chrom sweep" algorithm in BEDtools, making the assumption that 
// intervals are sorted by chrom and start.
//
// algorithm from https://github.com/arq5x/chrom_sweep/blob/master/chrom_sweep.py
//
// Assume that we will be operating on tibbles that are grouped by chromosome. 
// Eventually this wille enable efficient parallelization with multidplyr on the
// R side.
//

#include <Rcpp.h>
#include <queue>
using namespace Rcpp ;

// Only define fields for interval comparison. We will return a data frame
// that can be used to reconstruct any requested structure.
struct interval_t {
  int start;
  int end;
} Interval;

// store intervals under consideration
// std::queue <interval_t> IntervalCache ;

// empty an interval cache
void emptyCache( std::queue <interval_t> &q ) {
  std::queue<interval_t> empty ;
  std::swap( q, empty ) ;
}

// calculate overlap between two intervals
int intervalOverlap(interval_t a, interval_t b) {
    return std::min(a.end, b.end) - std::max(a.start, b.start);  
}

// is interval a after interval b?
bool intervalAfter(interval_t a, interval_t b) {
  if (a.start > b.end) {
    return true ;
  } else {
    return false ;
  } 
}

std::queue <interval_t> scanCache(interval_t curr_interval,
                                  std::queue <interval_t> interval_cache,
                                  std::queue <interval_t> intersection_cache) {

  std::queue <interval_t> temp_cache ;

  while ( !interval_cache.empty() ) {
    
    interval_t cached_interval = interval_cache.front() ;
    interval_cache.pop() ;
   
    if ( !intervalAfter(curr_interval, cached_interval) ) {
      
      temp_cache.push(curr_interval) ;
    
      if ( intervalOverlap(curr_interval, cached_interval) > 0 ) {
        intersection_cache.push(cached_interval) ;
      }
      
    } 
  } 
  
  return temp_cache ;      
}

std::queue <interval_t> createIntervals(Rcpp::DataFrame df) {

  std::queue <interval_t> intervals ;
  
  std::vector <int> df_start_v =
    Rcpp::as< std::vector <int> > (df["start"]); 
  
  std::vector <int> df_end_v =
    Rcpp::as< std::vector <int> > (df["end"]); 

  for (unsigned i = 0; i < df_start_v.size(); ++i) {
  
    interval_t interval;
    
    interval.start = df_start_v[i] ;
    interval.end = df_end_v[i] ;
     
    intervals.push(interval) ;
  }
  
  return intervals ;  
} 

void sweepIntervals(Rcpp::DataFrame df_a, Rcpp::DataFrame df_b) {
 
  // intersected intervals to output
  std::queue <interval_t> intersection_cache ;
  // intervals under consideration
  std::queue <interval_t> interval_cache ;
 
  std::queue <interval_t> df_a_intervals = createIntervals(df_a) ;
  std::queue <interval_t> df_b_intervals = createIntervals(df_b) ;

  interval_t curr_a_interval = df_a_intervals.front() ;
  interval_t curr_b_interval = df_b_intervals.front() ;
  
  df_a_intervals.pop() ;
  df_b_intervals.pop() ;
  
  while ( !df_a_intervals.empty() ) {
  
    interval_cache = scanCache(curr_b_interval, interval_cache, intersection_cache);
   
    while ( !df_b_intervals.empty() &&
            !intervalAfter(curr_a_interval, curr_b_interval) ) {
      
      if (intervalOverlap(curr_a_interval, curr_b_interval) > 0) {
        intersection_cache.push(curr_b_interval) ;  
      }
      
      interval_cache.push(curr_b_interval) ;
      
    }
    
    // XXX report the hits with the current a interval
    // reportHits(curr_a_interval, intersection_cache) ;
    
    // reset the intersection_cache
    emptyCache(intersection_cache) ;
    
    // get next b interval 
    interval_t curr_b_interval = df_b_intervals.front() ;
    df_b_intervals.pop() ;
  } 
}

// [[ Rcpp::export ]]

Rcpp::DataFrame intersect_(DataFrame df_a, DataFrame df_b) {
  
  // 1. sweep for intervals
  sweepIntervals(df_a, df_b) ;
    
  // 2. recreate intervals
  //std::vector <int> df_starts = intersection_cache.start
  //std::vector <int> df_ends = intersection_cache.end
  
   
//   return Rcpp::DataFrame::create( Named("start") = df_starts,
//                                   Named("end") = df_ends) ;
}
