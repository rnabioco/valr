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
#include <deque>
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
void emptyCache( std::queue <interval_t>& q ) {
  std::queue <interval_t> empty ;
  std::swap( q, empty ) ;
}

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

void scanCache(interval_t curr_interval,
               std::deque <interval_t> interval_cache,
               std::queue <interval_t>& intersection_cache) {

  while ( !interval_cache.empty() ) {
    
    interval_t cached_interval = interval_cache.front() ;
    interval_cache.pop_front() ;
   
    if ( !intervalAfter(curr_interval, cached_interval) ) {
    
      if ( intervalOverlap(curr_interval, cached_interval) > 0 ) {
        intersection_cache.push(cached_interval) ;
      }
      
    } 
  } 
}

std::queue <interval_t> createIntervals(Rcpp::DataFrame df) {

  std::vector <double> df_start_v =
    Rcpp::as< std::vector <double> >(df["start"]); 
  
  std::vector <double> df_end_v =
    Rcpp::as< std::vector <double> >(df["end"]); 

  std::queue <interval_t> intervals ;
 
  for (unsigned i = 0; i < df_start_v.size(); ++i) {
  
    interval_t interval;
    
    interval.start = df_start_v[i] ;
    interval.end = df_end_v[i] ;
     
    intervals.push(interval) ;
  }
  
  return intervals ;  
} 

void storeIntersections(interval_t query_interval,
                        std::queue <interval_t>& intersection_cache,
                        std::queue <intersection_t>& interval_intersections) {
  
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
}

// A = query
// B = database
std::queue <intersection_t> sweepIntervals(std::queue <interval_t> a_intervals,
                                           std::queue <interval_t> b_intervals) {
 
  // intersected intervals to output
  std::queue <interval_t> intersection_cache ;
  // intervals under consideration
  std::deque <interval_t> interval_cache ;
  // intersected intervals
  std::queue <intersection_t> interval_intersections ; 
 
  while ( !a_intervals.empty() ) {
  
    interval_t curr_a_interval = a_intervals.front() ;
    a_intervals.pop() ;
    
    // Rcpp::Rcout << "evaluating interval start: " << curr_a_interval.start << std::endl ;
    
    // this works?
    // scanCache(curr_a_interval, interval_cache, intersection_cache);
    Rcpp::Rcout << "cache size before scan: " << interval_cache.size() << std::endl ;
    scanCache(curr_a_interval, interval_cache, intersection_cache);
    Rcpp::Rcout << "cache size after scan: " << interval_cache.size() << std::endl ;
   
    while ( !b_intervals.empty() ) {
      
      interval_t curr_b_interval = b_intervals.front() ;
      b_intervals.pop() ;
  
      if (intervalAfter(curr_a_interval, curr_b_interval)) {
        break;
      }
      
      if (intervalOverlap(curr_a_interval, curr_b_interval) > 0) {
        intersection_cache.push(curr_b_interval) ;  
      }
      
      interval_cache.push_front(curr_b_interval) ;
      
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
    
    a_starts_v.push_back(intersection.a_start) ;
    a_ends_v.push_back(intersection.a_end) ;
    b_starts_v.push_back(intersection.b_start) ;
    b_ends_v.push_back(intersection.b_end) ;
    
  }
  
  return DataFrame::create( Named("a_start") = a_starts_v,
                            Named("a_end") = a_ends_v,
                            Named("b_start") = b_starts_v,
                            Named("b_end") = b_ends_v) ;
}
