//
// complement.cpp
//
// identify intervals not covered by input intervals
//

#include <Rcpp.h>
using namespace Rcpp ;

#include "intervals.h"
#include "genome.h"

interval_t
make_interval(std::string chrom, double start, double end) {
  
  interval_t interval ;
  
  interval.chrom = chrom ;
  interval.start = start ;
  interval.end = end ;
  
  return (interval) ;
}

std::list <interval_t>
complement_intervals(std::list<interval_t> intervals, std::map<std::string, double> genome) {

  std::list <interval_t> compl_intervals ;
  
  std::string last_chrom ; 
  double last_end ;
  
  int interval_count = 0;
  
  std::list<interval_t>::iterator it ; 
  for (it = intervals.begin(); it != intervals.end(); ++it) {
    
    interval_t curr_interval = *it ;
    
    // first interval on chrom
    if (interval_count == 0) {
 
      interval_t interval = make_interval(curr_interval.chrom, 1, curr_interval.start) ;   
      compl_intervals.push_back(interval) ; 
     
    } else if (curr_interval.chrom != last_chrom) {
      
      // switching chroms - add last interval on last chrom 
      double chrom_size = genome[last_chrom] ;
      interval_t last_interval = make_interval(last_chrom, last_end, chrom_size) ;
      compl_intervals.push_back(last_interval) ; 
      
      // add the first interval
      interval_t interval = make_interval(curr_interval.chrom, 1, curr_interval.start) ;   
      compl_intervals.push_back(interval) ; 
      
      interval_count = 0 ;
      
    } else {
      
      // internal interval on same chrom 
      interval_t interval = make_interval(curr_interval.chrom, last_end + 1, curr_interval.start) ;
      compl_intervals.push_back(interval) ; 
      
    }
    
    last_end = curr_interval.end ;
    last_chrom = curr_interval.chrom ;
    ++interval_count ;
 
  }  
  //  add final interval
  double chrom_size = genome[last_chrom] ;
  interval_t interval = make_interval(last_chrom, last_end, chrom_size) ;
  compl_intervals.push_back(interval) ;  
  
  return (compl_intervals) ;
}

// [[Rcpp::export]]
Rcpp::DataFrame
complement_impl(Rcpp::DataFrame interval_df, Rcpp::DataFrame genome_df) {
  
  std::list <interval_t> intervals = create_intervals(interval_df) ;
  std::map <std::string, double> genome = create_genome(genome_df) ;
  
  Rcpp::CharacterVector chroms_v ;
  Rcpp::NumericVector starts_v ;    
  Rcpp::NumericVector ends_v ;    
  
  std::list<interval_t> compl_intervals = complement_intervals(intervals, genome) ;

  std::list<interval_t>::const_iterator it; 
  for (it = compl_intervals.begin(); it != compl_intervals.end(); ++it) {
    
    interval_t interval = *it;
    
    chroms_v.push_back(interval.chrom) ; 
    starts_v.push_back(interval.start) ;
    ends_v.push_back(interval.end) ;
    
  }
  
  return DataFrame::create( Named("chrom") = chroms_v,
                            Named("start") = starts_v,
                            Named("end") = ends_v) ;
      
}
