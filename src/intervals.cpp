#include <Rcpp.h>
#include "intervals.h"

// calculate overlap between two intervals
int
interval_overlap(const interval_t a, const interval_t b) {
  return std::min(a.end, b.end) - std::max(a.start, b.start);  
}

// is interval a after interval b?
bool
interval_after(const interval_t a, const interval_t b) {
  if (a.start > b.end) {
    return true ;
  } else {
    return false ;
  } 
}

std::list <interval_t>
create_intervals(Rcpp::DataFrame df) {
 
  std::vector <std::string> df_chrom_v = 
    Rcpp::as< std::vector <std::string> >(df["chrom"]); 
  
  std::vector <int> df_start_v =
    Rcpp::as< std::vector <int> >(df["start"]); 
  
  std::vector <int> df_end_v =
    Rcpp::as< std::vector <int> >(df["end"]); 
  
  std::list <interval_t> intervals ;
  
  for (unsigned i = 0; i < df_start_v.size(); ++i) {
    
    interval_t interval;
    
    interval.chrom = df_chrom_v[i] ;   
    interval.start = df_start_v[i] ;
    interval.end = df_end_v[i] ;
    
    intervals.push_back(interval) ;
  }
  
  return intervals ;  
} 
