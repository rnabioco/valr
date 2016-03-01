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
save_interval(std::list<interval_t>& intervals, std::string chrom, double start, double end) {
  interval_t interval = make_interval(chrom, start, end) ;
  intervals.push_back(interval) ;
}

std::list <interval_t>
complement_intervals(std::list<interval_t> intervals, std::map<std::string, double> genome) {

  std::list <interval_t> compl_intervals ;
  
  std::string prev_chrom = ""; 
  double prev_end = 0;
  
  int interval_count = 0;
  
  std::list<interval_t>::iterator it ; 
  for (it = intervals.begin(); it != intervals.end(); ++it) {
    
    interval_t curr_interval = *it ;
    
    // first interval on chrom
    if (interval_count == 0) {

      if (curr_interval.start > 1) { 
        save_interval(compl_intervals, curr_interval.chrom, 1, curr_interval.start) ;
      }
      
    } else if (prev_chrom != "" && curr_interval.chrom != prev_chrom) {
      
      // switching chroms - add last interval on previous chrom 
      double chrom_size = genome[prev_chrom] ;
      save_interval(compl_intervals, prev_chrom, prev_end + 1, chrom_size) ;
      
      // add the first interval
      if (curr_interval.start > 1) { 
        save_interval(compl_intervals, curr_interval.chrom, 1, curr_interval.start) ;
      } 
      
      interval_count = 0 ;
      
    } else {
      // internal interval on same chrom 
      save_interval(compl_intervals, curr_interval.chrom, prev_end + 1, curr_interval.start) ;
    }
    
    prev_end = curr_interval.end ;
    prev_chrom = curr_interval.chrom ;
    ++interval_count ;
 
  }  
  
  //  add final interval
  double chrom_size = genome[prev_chrom] ;
  if (prev_end < chrom_size) {
    save_interval(compl_intervals, prev_chrom, prev_end + 1, chrom_size) ;
  } 
  
  return compl_intervals ;
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
