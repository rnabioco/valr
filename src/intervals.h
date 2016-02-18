#ifndef INTERVALS_INCLUDE

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

#define INTERVALS_INCLUDE
#endif /* INTERVALS_INCLUDE */