// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/valr.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bed12toexons_impl
DataFrame bed12toexons_impl(DataFrame x);
RcppExport SEXP _valr_bed12toexons_impl(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bed12toexons_impl(x));
    return rcpp_result_gen;
END_RCPP
}
// closest_impl
DataFrame closest_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y, IntegerVector grp_idx_x, IntegerVector grp_idx_y, const std::string& suffix_x, const std::string& suffix_y);
RcppExport SEXP _valr_closest_impl(SEXP xSEXP, SEXP ySEXP, SEXP grp_idx_xSEXP, SEXP grp_idx_ySEXP, SEXP suffix_xSEXP, SEXP suffix_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grp_idx_x(grp_idx_xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grp_idx_y(grp_idx_ySEXP);
    Rcpp::traits::input_parameter< const std::string& >::type suffix_x(suffix_xSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type suffix_y(suffix_ySEXP);
    rcpp_result_gen = Rcpp::wrap(closest_impl(x, y, grp_idx_x, grp_idx_y, suffix_x, suffix_y));
    return rcpp_result_gen;
END_RCPP
}
// complement_impl
DataFrame complement_impl(ValrGroupedDataFrame gdf, DataFrame genome);
RcppExport SEXP _valr_complement_impl(SEXP gdfSEXP, SEXP genomeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type gdf(gdfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type genome(genomeSEXP);
    rcpp_result_gen = Rcpp::wrap(complement_impl(gdf, genome));
    return rcpp_result_gen;
END_RCPP
}
// coverage_impl
DataFrame coverage_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y, IntegerVector x_grp_indexes, IntegerVector y_grp_indexes);
RcppExport SEXP _valr_coverage_impl(SEXP xSEXP, SEXP ySEXP, SEXP x_grp_indexesSEXP, SEXP y_grp_indexesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x_grp_indexes(x_grp_indexesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y_grp_indexes(y_grp_indexesSEXP);
    rcpp_result_gen = Rcpp::wrap(coverage_impl(x, y, x_grp_indexes, y_grp_indexes));
    return rcpp_result_gen;
END_RCPP
}
// dist_impl
DataFrame dist_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y, IntegerVector x_grp_indexes, IntegerVector y_grp_indexes, std::string distcalc);
RcppExport SEXP _valr_dist_impl(SEXP xSEXP, SEXP ySEXP, SEXP x_grp_indexesSEXP, SEXP y_grp_indexesSEXP, SEXP distcalcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x_grp_indexes(x_grp_indexesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y_grp_indexes(y_grp_indexesSEXP);
    Rcpp::traits::input_parameter< std::string >::type distcalc(distcalcSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_impl(x, y, x_grp_indexes, y_grp_indexes, distcalc));
    return rcpp_result_gen;
END_RCPP
}
// gcoverage_impl
DataFrame gcoverage_impl(const ValrGroupedDataFrame& gdf, const IntegerVector& max_coords);
RcppExport SEXP _valr_gcoverage_impl(SEXP gdfSEXP, SEXP max_coordsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ValrGroupedDataFrame& >::type gdf(gdfSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type max_coords(max_coordsSEXP);
    rcpp_result_gen = Rcpp::wrap(gcoverage_impl(gdf, max_coords));
    return rcpp_result_gen;
END_RCPP
}
// intersect_impl
DataFrame intersect_impl(ValrGroupedDataFrame x, ValrGroupedDataFrame y, IntegerVector x_grp_indexes, IntegerVector y_grp_indexes, bool invert, const std::string& suffix_x, const std::string& suffix_y);
RcppExport SEXP _valr_intersect_impl(SEXP xSEXP, SEXP ySEXP, SEXP x_grp_indexesSEXP, SEXP y_grp_indexesSEXP, SEXP invertSEXP, SEXP suffix_xSEXP, SEXP suffix_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x_grp_indexes(x_grp_indexesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y_grp_indexes(y_grp_indexesSEXP);
    Rcpp::traits::input_parameter< bool >::type invert(invertSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type suffix_x(suffix_xSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type suffix_y(suffix_ySEXP);
    rcpp_result_gen = Rcpp::wrap(intersect_impl(x, y, x_grp_indexes, y_grp_indexes, invert, suffix_x, suffix_y));
    return rcpp_result_gen;
END_RCPP
}
// merge_impl
DataFrame merge_impl(ValrGroupedDataFrame gdf, int max_dist, bool collapse);
RcppExport SEXP _valr_merge_impl(SEXP gdfSEXP, SEXP max_distSEXP, SEXP collapseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type gdf(gdfSEXP);
    Rcpp::traits::input_parameter< int >::type max_dist(max_distSEXP);
    Rcpp::traits::input_parameter< bool >::type collapse(collapseSEXP);
    rcpp_result_gen = Rcpp::wrap(merge_impl(gdf, max_dist, collapse));
    return rcpp_result_gen;
END_RCPP
}
// partition_impl
DataFrame partition_impl(const ValrGroupedDataFrame& gdf, int max_dist);
RcppExport SEXP _valr_partition_impl(SEXP gdfSEXP, SEXP max_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const ValrGroupedDataFrame& >::type gdf(gdfSEXP);
    Rcpp::traits::input_parameter< int >::type max_dist(max_distSEXP);
    rcpp_result_gen = Rcpp::wrap(partition_impl(gdf, max_dist));
    return rcpp_result_gen;
END_RCPP
}
// shuffle_impl
DataFrame shuffle_impl(DataFrame df, DataFrame incl, bool within, int max_tries, int seed);
RcppExport SEXP _valr_shuffle_impl(SEXP dfSEXP, SEXP inclSEXP, SEXP withinSEXP, SEXP max_triesSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type incl(inclSEXP);
    Rcpp::traits::input_parameter< bool >::type within(withinSEXP);
    Rcpp::traits::input_parameter< int >::type max_tries(max_triesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(shuffle_impl(df, incl, within, max_tries, seed));
    return rcpp_result_gen;
END_RCPP
}
// subtract_impl
DataFrame subtract_impl(ValrGroupedDataFrame gdf_x, ValrGroupedDataFrame gdf_y, IntegerVector x_grp_indexes, IntegerVector y_grp_indexes);
RcppExport SEXP _valr_subtract_impl(SEXP gdf_xSEXP, SEXP gdf_ySEXP, SEXP x_grp_indexesSEXP, SEXP y_grp_indexesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type gdf_x(gdf_xSEXP);
    Rcpp::traits::input_parameter< ValrGroupedDataFrame >::type gdf_y(gdf_ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x_grp_indexes(x_grp_indexesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y_grp_indexes(y_grp_indexesSEXP);
    rcpp_result_gen = Rcpp::wrap(subtract_impl(gdf_x, gdf_y, x_grp_indexes, y_grp_indexes));
    return rcpp_result_gen;
END_RCPP
}
