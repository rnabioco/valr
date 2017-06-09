#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP valr_absdist_impl(SEXP, SEXP);
extern SEXP valr_bed12toexons_impl(SEXP);
extern SEXP valr_closest_impl(SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_complement_impl(SEXP, SEXP);
extern SEXP valr_coverage_impl(SEXP, SEXP);
extern SEXP valr_flank_impl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_intersect_impl(SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_makewindows_impl(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_merge_impl(SEXP, SEXP, SEXP);
extern SEXP valr_random_impl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_reldist_impl(SEXP, SEXP);
extern SEXP valr_shuffle_impl(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP valr_subtract_impl(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"valr_absdist_impl",      (DL_FUNC) &valr_absdist_impl,      2},
  {"valr_bed12toexons_impl", (DL_FUNC) &valr_bed12toexons_impl, 1},
  {"valr_closest_impl",      (DL_FUNC) &valr_closest_impl,      4},
  {"valr_complement_impl",   (DL_FUNC) &valr_complement_impl,   2},
  {"valr_coverage_impl",     (DL_FUNC) &valr_coverage_impl,     2},
  {"valr_flank_impl",        (DL_FUNC) &valr_flank_impl,        8},
  {"valr_intersect_impl",    (DL_FUNC) &valr_intersect_impl,    4},
  {"valr_makewindows_impl",  (DL_FUNC) &valr_makewindows_impl,  5},
  {"valr_merge_impl",        (DL_FUNC) &valr_merge_impl,        3},
  {"valr_random_impl",       (DL_FUNC) &valr_random_impl,       6},
  {"valr_reldist_impl",      (DL_FUNC) &valr_reldist_impl,      2},
  {"valr_shuffle_impl",      (DL_FUNC) &valr_shuffle_impl,      5},
  {"valr_subtract_impl",     (DL_FUNC) &valr_subtract_impl,     2},
  {NULL, NULL, 0}
};

void R_init_valr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
