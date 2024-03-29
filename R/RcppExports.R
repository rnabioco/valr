# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

bed12toexons_impl <- function(x) {
    .Call(`_valr_bed12toexons_impl`, x)
}

closest_impl <- function(x, y, grp_idx_x, grp_idx_y, suffix_x, suffix_y) {
    .Call(`_valr_closest_impl`, x, y, grp_idx_x, grp_idx_y, suffix_x, suffix_y)
}

complement_impl <- function(gdf, genome) {
    .Call(`_valr_complement_impl`, gdf, genome)
}

coverage_impl <- function(x, y, x_grp_indexes, y_grp_indexes) {
    .Call(`_valr_coverage_impl`, x, y, x_grp_indexes, y_grp_indexes)
}

dist_impl <- function(x, y, x_grp_indexes, y_grp_indexes, distcalc) {
    .Call(`_valr_dist_impl`, x, y, x_grp_indexes, y_grp_indexes, distcalc)
}

flank_impl <- function(df, genome, both = 0, left = 0, right = 0, fraction = FALSE, stranded = FALSE, trim = FALSE) {
    .Call(`_valr_flank_impl`, df, genome, both, left, right, fraction, stranded, trim)
}

gcoverage_impl <- function(gdf, max_coords) {
    .Call(`_valr_gcoverage_impl`, gdf, max_coords)
}

intersect_impl <- function(x, y, x_grp_indexes, y_grp_indexes, invert = FALSE, suffix_x = ".x", suffix_y = ".y") {
    .Call(`_valr_intersect_impl`, x, y, x_grp_indexes, y_grp_indexes, invert, suffix_x, suffix_y)
}

makewindows_impl <- function(df, win_size = 0L, num_win = 0L, step_size = 0L, reverse = FALSE) {
    .Call(`_valr_makewindows_impl`, df, win_size, num_win, step_size, reverse)
}

merge_impl <- function(gdf, max_dist = 0L, collapse = TRUE) {
    .Call(`_valr_merge_impl`, gdf, max_dist, collapse)
}

partition_impl <- function(gdf, max_dist = -1L) {
    .Call(`_valr_partition_impl`, gdf, max_dist)
}

random_impl <- function(genome, length, n, seed = 0L) {
    .Call(`_valr_random_impl`, genome, length, n, seed)
}

shuffle_impl <- function(df, incl, within = FALSE, max_tries = 1000L, seed = 0L) {
    .Call(`_valr_shuffle_impl`, df, incl, within, max_tries, seed)
}

subtract_impl <- function(gdf_x, gdf_y, x_grp_indexes, y_grp_indexes) {
    .Call(`_valr_subtract_impl`, gdf_x, gdf_y, x_grp_indexes, y_grp_indexes)
}

