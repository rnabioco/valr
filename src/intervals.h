#ifndef INTERVALS_INCLUDE
#define INTERVALS_INCLUDE

struct interval_t {
  std::string chrom ;
  int start, end ;
};

struct intersection_t {
  interval_t a, b ;
};

// calculate overlap between two intervals
int intervalOverlap(const interval_t a, const interval_t b) ;

// is interval a after interval b?
bool intervalAfter(const interval_t a, const interval_t b) ;

std::list <interval_t> createIntervals(Rcpp::DataFrame df) ;

#endif /* INTERVALS_INCLUDE */