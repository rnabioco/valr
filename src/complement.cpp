//
// complement.cpp
//
// identify intervals not covered by input intervals
//

#include <Rcpp.h>
using namespace Rcpp ;

#include "intervals.h"
#include "genome.h"

void
save_interval(std::list<interval_t>& intervals, std::string chrom, int start, int end) {
  interval_t interval = make_interval(chrom, start, end) ;
  intervals.push_back(interval) ;
}

std::list <interval_t>
complement_intervals(std::list<interval_t> intervals, std::map<std::string, int> genome) {

  interval_t prev_interval = make_interval("", 0, 0) ;
  std::list <interval_t> compl_intervals ;
  
  std::list<interval_t>::iterator it ; 
  for (it = intervals.begin(); it != intervals.end(); ++it) {
    
    interval_t curr_interval = *it ;
    
    // first interval on chrom
    if ( it == intervals.begin() ) {

      if (curr_interval.start > 1) { 
        save_interval(compl_intervals, curr_interval.chrom, 1, curr_interval.start) ;
      }
      
    } else if (prev_interval.chrom != "" && curr_interval.chrom != prev_interval.chrom) {
      
      // switching chroms - add last interval on previous chrom 
      int chrom_size = genome[prev_interval.chrom] ;
      save_interval(compl_intervals, prev_interval.chrom, prev_interval.end + 1, chrom_size) ;
      
      // add the first interval
      if (curr_interval.start > 1) { 
        save_interval(compl_intervals, curr_interval.chrom, 1, curr_interval.start) ;
      } 
      
    } else if (curr_interval.chrom == prev_interval.chrom ) {
      
      // internal interval on same chrom 
      save_interval(compl_intervals, curr_interval.chrom, prev_interval.end + 1, curr_interval.start) ;
      
    } else if (it == intervals.end() ) {
      
       //  add final interval
      int chrom_size = genome[curr_interval.chrom] ;
      if (curr_interval.end < chrom_size) {
        save_interval(compl_intervals, curr_interval.chrom, curr_interval.end, chrom_size) ;
      } 
    }
    
    prev_interval = *it ;
  }  
 
  return compl_intervals ;
}

// [[Rcpp::export]]
Rcpp::DataFrame
complement_impl(Rcpp::DataFrame interval_df, Rcpp::DataFrame genome_df) {
  
  std::list <interval_t> intervals = create_intervals(interval_df) ;
  std::map <std::string, int> genome = create_genome(genome_df) ;
  
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
