// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inner_project_loop
NumericMatrix inner_project_loop(int no_sp, int no_w, NumericMatrix n, NumericMatrix A, NumericMatrix B, NumericMatrix S, NumericVector w_min_idx);
RcppExport SEXP _mizer_inner_project_loop(SEXP no_spSEXP, SEXP no_wSEXP, SEXP nSEXP, SEXP ASEXP, SEXP BSEXP, SEXP SSEXP, SEXP w_min_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_sp(no_spSEXP);
    Rcpp::traits::input_parameter< int >::type no_w(no_wSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_min_idx(w_min_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_project_loop(no_sp, no_w, n, A, B, S, w_min_idx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mizer_inner_project_loop", (DL_FUNC) &_mizer_inner_project_loop, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_mizer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
