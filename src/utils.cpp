#include "valr.h"


int intervalOverlap(const interval_t& a, const interval_t& b) {
  return(std::min(a.stop, b.stop) - std::max(a.start, b.start)) ;    
}

// the value field of intervals in the returned vector correspond to the index
// of the interval in the original dataframe (i.e., the values of the
// SlicingIndex)
intervalVector makeIntervalVector(DataFrame df, SlicingIndex si) {
  
  intervalVector intervals ;
  
  IntegerVector starts = df["start"] ;
  IntegerVector ends   = df["end"] ;
  
  int size = si.size() ;
  
  for( int i=0; i<size; ++i) {
    int j = si[i] ;
    intervals.push_back(interval_t(starts[j], ends[j], j)) ;
  }
  return intervals ;
}

icl_interval_set_t makeIclIntervalSet(DataFrame df, SlicingIndex indices) {

  icl_interval_set_t iset ;
  
  IntegerVector starts = df["start"] ;
  IntegerVector ends   = df["end"] ;
  
  int size = indices.size() ;
  
  for( int i=0; i<size; ++i) {
    int j = indices[i] ;
    icl_interval_t intvl = icl_interval_t::closed(starts[j], ends[j]) ;
    iset.add(intvl) ;
  }

  return iset ;
}
 
 // generic function to loop through x and y bed_tbls and apply a function to each matched chromosomes
 
// void chromLoop( const GroupedDataFrame& x, const GroupedDataFrame& y,
//                std::function<void (intervalVector&, intervalVector&)> fxn) {

//   auto data_x = x.data() ;
//   auto data_y = y.data() ;
   
   //   auto ng_x = x.ngroups() ;
   //   auto ng_y = y.ngroups() ;
   
   //CharacterVector chrom_x = data_x["chrom"];
   //CharacterVector chrom_y = data_y["chrom"];
   
   
   //   GroupedDataFrame::group_iterator git_x = x.group_begin() ;
   // for( int nx=0; nx<ng_x; nx++, ++git_x){
     
     // SlicingIndex gi_x = *git_x ;
     //   int group_x = gi_x[0];
     
     // GroupedDataFrame::group_iterator git_y = y.group_begin() ;
     //     for( int ny=0; ny<ng_y; ny++, ++git_y) {
       
       //SlicingIndex gi_y = *git_y ;
       // int group_y = gi_y[0];
       
       // if( chrom_x[group_x] == chrom_y[group_y] ) {
         
         //intervalVector vx = makeIntervalVector(data_x, gi_x) ;
         // intervalVector vy = makeIntervalVector(data_y, gi_y) ;
         
         //fxn(vx, vy ) ; 
         //}
       //} 
     //}
   // }

