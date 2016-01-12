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
} Interval ;

struct intersection_t {
  int a_start;
  int a_end;
  int b_start;
  int b_end;
} Intersection ;

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

std::queue <interval_t> createIntervals(DataFrame df) {

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

std::queue <intersection_t> storeIntersections(interval_t query_interval,
                                               std::queue <interval_t> intersection_cache,
                                               std::queue <intersection_t> interval_intersections) {
  while ( !intersection_cache.empty() ) {
    
    interval_t cache_interval = intersection_cache.front() ;
    intersection_cache.pop() ;
  
    // create the intersection 
    intersection_t intersection ;
    
    intersection.a_start = query_interval.start ;
    intersection.a_end = query_interval.end ;
    intersection.b_start = cache_interval.start ;
    intersection.b_end = cache_interval.end ;
    
    interval_intersections.push(intersection) ;
  }
  
  return interval_intersections ;
}

std::queue <intersection_t> sweepIntervals(std::queue <interval_t> a_intervals,
                                           std::queue <interval_t> b_intervals) {
 
  // intersected intervals to output
  std::queue <interval_t> intersection_cache ;
  // intervals under consideration
  std::queue <interval_t> interval_cache ;
  // intersected intervals
  std::queue <intersection_t> interval_intersections ; 
 
  interval_t curr_a_interval = a_intervals.front() ;
  interval_t curr_b_interval = b_intervals.front() ;
  
  a_intervals.pop() ;
  b_intervals.pop() ;
  
  while ( !a_intervals.empty() ) {
  
    interval_cache = scanCache(curr_b_interval, interval_cache, intersection_cache);
   
    while ( !b_intervals.empty() &&
            !intervalAfter(curr_a_interval, curr_b_interval) ) {
      
      if (intervalOverlap(curr_a_interval, curr_b_interval) > 0) {
        intersection_cache.push(curr_b_interval) ;  
      }
      
      interval_cache.push(curr_b_interval) ;
      
    }
    
    // store hits with the current a interval
    interval_intersections = storeIntersections(curr_a_interval,
                                                intersection_cache,
                                                interval_intersections) ;
    
    // reset the intersection_cache - is it already empty?
    emptyCache(intersection_cache) ;
    
    // get next b interval 
    interval_t curr_b_interval = b_intervals.front() ;
    b_intervals.pop() ;
  } 
  
  return interval_intersections ;
}

// [[Rcpp::export]]
DataFrame intersect_cpp(DataFrame df_a, DataFrame df_b) {
  
  std::queue <interval_t> a_intervals = createIntervals(df_a) ;
  std::queue <interval_t> b_intervals = createIntervals(df_b) ;
 
  std::queue <intersection_t> interval_intersections =
    sweepIntervals(a_intervals, b_intervals) ;

  Rcpp::NumericVector a_starts_v ;    
  Rcpp::NumericVector a_ends_v ;    
  Rcpp::NumericVector b_starts_v ;    
  Rcpp::NumericVector b_ends_v ;    
 
  while ( !interval_intersections.empty()) {
    
    intersection_t intersection = interval_intersections.front() ;
    interval_intersections.pop() ;
    
    a_starts_v.push_front(intersection.a_start) ;
    a_ends_v.push_front(intersection.a_end) ;
    b_starts_v.push_front(intersection.b_start) ;
    b_ends_v.push_front(intersection.b_end) ;
    
  }
  
  return DataFrame::create( Named("a_start") = a_starts_v,
                            Named("a_end") = a_ends_v,
                            Named("b_start") = b_starts_v,
                            Named("b_end") = b_ends_v) ;
}
