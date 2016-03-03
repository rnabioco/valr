#ifndef Rbedtools_intervals_H
#define Rbedtools_intervals_H 

struct interval_t {
  std::string chrom ;
  int start, end ;
};

struct intersection_t {
  interval_t a, b ;
};

// create a new interval
interval_t make_interval(std::string chrom, int start, int end) ;
    
// calculate overlap between two intervals
int interval_overlap(const interval_t a, const interval_t b) ;

// is interval a after interval b?
bool interval_after(const interval_t a, const interval_t b) ;

std::list <interval_t> create_intervals(Rcpp::DataFrame df) ;

#endif 
