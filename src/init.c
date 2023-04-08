#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _valr_bed12toexons_impl(void *);
extern SEXP _valr_closest_impl(void *, void *, void *, void *, void *, void *);
extern SEXP _valr_complement_impl(void *, void *);
extern SEXP _valr_coverage_impl(void *, void *, void *, void *);
extern SEXP _valr_dist_impl(void *, void *, void *, void *, void *);
extern SEXP _valr_flank_impl(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _valr_intersect_impl(void *, void *, void *, void *, void *, void *, void *);
extern SEXP _valr_makewindows_impl(void *, void *, void *, void *, void *);
extern SEXP _valr_merge_impl(void *, void *, void *);
extern SEXP _valr_partition_impl(void *, void *);
extern SEXP _valr_random_impl(void *, void *, void *, void *);
extern SEXP _valr_shuffle_impl(void *, void *, void *, void *, void *);
extern SEXP _valr_subtract_impl(void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_valr_bed12toexons_impl", (DL_FUNC) &_valr_bed12toexons_impl, 1},
    {"_valr_closest_impl",      (DL_FUNC) &_valr_closest_impl,      6},
    {"_valr_complement_impl",   (DL_FUNC) &_valr_complement_impl,   2},
    {"_valr_coverage_impl",     (DL_FUNC) &_valr_coverage_impl,     4},
    {"_valr_dist_impl",         (DL_FUNC) &_valr_dist_impl,         5},
    {"_valr_flank_impl",        (DL_FUNC) &_valr_flank_impl,        8},
    {"_valr_intersect_impl",    (DL_FUNC) &_valr_intersect_impl,    7},
    {"_valr_makewindows_impl",  (DL_FUNC) &_valr_makewindows_impl,  5},
    {"_valr_merge_impl",        (DL_FUNC) &_valr_merge_impl,        3},
    {"_valr_partition_impl",    (DL_FUNC) &_valr_partition_impl,    2},
    {"_valr_random_impl",       (DL_FUNC) &_valr_random_impl,       4},
    {"_valr_shuffle_impl",      (DL_FUNC) &_valr_shuffle_impl,      5},
    {"_valr_subtract_impl",     (DL_FUNC) &_valr_subtract_impl,     4},
    {NULL, NULL, 0}
};

void R_init_valr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
