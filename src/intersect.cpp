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
  char chrom ;
  int start, end ;
} Interval ;

struct intersection_t {
  interval_t a, b ;
} Intersection ;

// calculate overlap between two intervals
int intervalOverlap(const interval_t a, const interval_t b) {
  return std::min(a.end, b.end) - std::max(a.start, b.start);  
}

// is interval a after interval b?
bool intervalAfter(const interval_t a, const interval_t b) {
  if (a.start > b.end) {
    return true ;
  } else {
    return false ;
  } 
}


std::list <interval_t> createIntervals(Rcpp::DataFrame df) {
  
  std::vector <double> df_start_v =
    Rcpp::as< std::vector <double> >(df["start"]); 
  
  std::vector <double> df_end_v =
    Rcpp::as< std::vector <double> >(df["end"]); 
  
  std::list <interval_t> intervals ;
  
  for (unsigned i = 0; i < df_start_v.size(); ++i) {
    
    interval_t interval;
    
    interval.start = df_start_v[i] ;
    interval.end = df_end_v[i] ;
    
    intervals.push_back(interval) ;
  }
  
  return intervals ;  
} 

void storeIntersections(interval_t query_interval,
                        std::queue <interval_t>& intersection_cache,
                        std::list <intersection_t>& interval_intersections) {
  
  while ( !intersection_cache.empty() ) {
    
    interval_t cache_interval = intersection_cache.front() ;
    intersection_cache.pop() ;
    
    // create the intersection 
    intersection_t intersection ;
    
    intersection.a.start = query_interval.start ;
    intersection.a.end = query_interval.end ;
    intersection.b.start = cache_interval.start ;
    intersection.b.end = cache_interval.end ;
    
    interval_intersections.push_back(intersection) ;
  }
}

//
std::list <interval_t> scanCache(interval_t curr_interval,
                                 std::list <interval_t> interval_cache,
                                 std::queue <interval_t>& intersection_cache) {
  
  std::list <interval_t> temp_cache ;
  
  std::list<interval_t>::const_iterator it; 
  for (it = interval_cache.begin(); it != interval_cache.end(); ++it) {
    
    interval_t cached_interval = *it ;
    
    if ( !intervalAfter(curr_interval, cached_interval) ) {
      
      temp_cache.push_back(cached_interval) ;
      
      if ( intervalOverlap(curr_interval, cached_interval) > 0 ) {
        intersection_cache.push(cached_interval) ;
      }
      
    } 
  } 
  return temp_cache ;
}

// A = query
// B = database
std::list <intersection_t> sweepIntervals(std::list <interval_t> a_intervals,
                                          std::list <interval_t> b_intervals) {
  
  // cached intersections
  std::queue <interval_t> intersection_cache ;
  // intervals under consideration
  std::list <interval_t> interval_cache ;
  // intersected intervals
  std::list <intersection_t> interval_intersections ; 
  
  std::list<interval_t>::iterator a_it ; 
  for (a_it = a_intervals.begin(); a_it != a_intervals.end(); ++a_it) {
    
    interval_t curr_a_interval = *a_it;
    
    interval_cache = scanCache(curr_a_interval, interval_cache, intersection_cache);
    
    while ( !b_intervals.empty() ) {
      
      interval_t curr_b_interval = b_intervals.front() ;
      
      if ( !intervalAfter(curr_a_interval, curr_b_interval) ) { 
        
        if (intervalOverlap(curr_a_interval, curr_b_interval) > 0) {
          intersection_cache.push(curr_b_interval) ;  
        }
        
        interval_cache.push_front(curr_b_interval) ;
        
        b_intervals.pop_front() ;
      } 
    }
    
    // store hits with the current a interval
    storeIntersections(curr_a_interval,
                       intersection_cache,
                       interval_intersections) ;
    
  } 
  
  return interval_intersections ;
}

// [[Rcpp::export]]
DataFrame intersect_cpp(DataFrame df_a, DataFrame df_b) {
  
  std::list <interval_t> a_intervals = createIntervals(df_a) ;
  std::list <interval_t> b_intervals = createIntervals(df_b) ;
  
  std::list <intersection_t> interval_intersections =
    sweepIntervals(a_intervals, b_intervals) ;
  
  Rcpp::NumericVector a_starts_v ;    
  Rcpp::NumericVector a_ends_v ;    
  Rcpp::NumericVector b_starts_v ;    
  Rcpp::NumericVector b_ends_v ;    
  
  std::list<intersection_t>::const_iterator it; 
  for (it = interval_intersections.begin(); it != interval_intersections.end(); ++it) {
    
    intersection_t intersection = *it;
    
    a_starts_v.push_back(intersection.a.start) ;
    a_ends_v.push_back(intersection.a.end) ;
    b_starts_v.push_back(intersection.b.start) ;
    b_ends_v.push_back(intersection.b.end) ;
    
  }
  
  return DataFrame::create( Named("a_start") = a_starts_v,
                            Named("a_end") = a_ends_v,
                            Named("b_start") = b_starts_v,
                            Named("b_end") = b_ends_v) ;
}
