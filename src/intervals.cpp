#include <Rcpp.h>
#include "Rbedtools.h"

// calculate overlap between two intervals
// If >0, the number of bp of overlap
// If 0,  they are book-ended.
// If <0, the distance in bp between them
// 
//[[Rcpp::export]]
int interval_overlap(int const& start_x, int const& end_x, int const& start_y, int const& end_y) {
  return(std::min(end_x, end_y) - std::max(start_x, start_y)) ;    
}


