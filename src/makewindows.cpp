#include "valr.h"

IntegerVector rcppSeq(int from, int to, int by = 1) {
  //determine output size
  std::size_t n = ((to - from) / by) + 1 ;
  //generate output vector
  IntegerVector res = Rcpp::rep(from, n) ;
  // iterate through from to to by step size
  int idx = 0 ;
  for (int i = from; i <= to; i += by ){
    res[idx] = i ;
    ++idx ;
  }
  return res ;
}

//[[Rcpp::export]]
DataFrame makewindows_impl(DataFrame df,
                           int step_size = 0,
                           bool reverse = false,
                           std::string col_start = "start",
                           std::string col_end = "end",
                           std::string col_win_size = ".win_size") {

  NumericVector starts = df[col_start] ;
  NumericVector ends = df[col_end] ;
  NumericVector win_sizes = df[col_win_size] ;

  // all vectors are the same size, take size from starts
  size_t n = starts.size() ;

  std::vector<int> starts_out ;
  std::vector<int> ends_out ;
  std::vector<int> df_idxs ;
  std::vector<int> window_ids;

  for(int i = 0; i < n ; ++i){
    auto start = starts[i] ;
    auto end = ends[i] ;
    auto win_size = win_sizes[i] ;
    int by = win_size - step_size ;

    // iterate by step_size
    auto ivl_starts_out = rcppSeq(start, end, by) ;

    // allocate new ends bases on new starts
    int out_size = ivl_starts_out.size() ;

    IntegerVector x_idxs = Rcpp::rep(i, out_size) ;
    IntegerVector ivl_ends_out(out_size) ;

    // make window ids
    std::vector<int> ids ;

    for(int i = 0; i < out_size ; ++i){
      // make end equal to end if past end
      ivl_ends_out[i] = ivl_starts_out[i] + win_size < end ? \
                        ivl_starts_out[i] + win_size : end ;
      ids.push_back(i + 1) ;
    }

    // remove  ivls if start == end
    for(int i = 0; i < out_size ; ++i){
      if(ivl_starts_out[i] == ivl_ends_out[i]){
        x_idxs.erase(x_idxs.begin() + i) ;
        ids.erase(ids.begin() + i) ;
        ivl_starts_out.erase(ivl_starts_out.begin() + i) ;
        ivl_ends_out.erase(ivl_ends_out.begin() + i) ;
      }
    }
    if(reverse){
      std::reverse(ids.begin(), ids.end());
    }

    // add new ivls
    starts_out.insert(starts_out.end(), ivl_starts_out.begin(), ivl_starts_out.end());
    ends_out.insert(ends_out.end(), ivl_ends_out.begin(), ivl_ends_out.end());
    df_idxs.insert(df_idxs.end(), x_idxs.begin(), x_idxs.end());
    window_ids.insert(window_ids.end(), ids.begin(), ids.end());
  }

  DataFrame out = DataFrameSubsetVisitors(df, names(df)).subset(df_idxs, "data.frame");
  out[col_start] = starts_out ;
  out[col_end] = ends_out ;
  out[".win_id"] = window_ids ;

  return out ;
}

/*** R
x <- trbl_interval(
  ~chrom, ~start, ~end, ~name, ~score, ~strand, ~.row_id, ~.win_size,
  "chr1", 100,    200,  'A',   '.',    '+', 1, 10
)

makewindows_impl(x) %>% as_data_frame()
*/
