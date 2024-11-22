// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calc_HOT
List calc_HOT(IntegerMatrix G, int kid, IntegerVector par);
RcppExport SEXP _MixedUpParents_calc_HOT(SEXP GSEXP, SEXP kidSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type kid(kidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_HOT(G, kid, par));
    return rcpp_result_gen;
END_RCPP
}
// intersect_ancestry_intervals_rcpp
IntegerMatrix intersect_ancestry_intervals_rcpp(IntegerVector grp1, IntegerVector grp2, int nC, IntegerMatrix V1, IntegerMatrix V2, IntegerMatrix X1, IntegerMatrix X2, IntegerVector nv1, IntegerVector nv2);
RcppExport SEXP _MixedUpParents_intersect_ancestry_intervals_rcpp(SEXP grp1SEXP, SEXP grp2SEXP, SEXP nCSEXP, SEXP V1SEXP, SEXP V2SEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP nv1SEXP, SEXP nv2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type grp1(grp1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type grp2(grp2SEXP);
    Rcpp::traits::input_parameter< int >::type nC(nCSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nv1(nv1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nv2(nv2SEXP);
    rcpp_result_gen = Rcpp::wrap(intersect_ancestry_intervals_rcpp(grp1, grp2, nC, V1, V2, X1, X2, nv1, nv2));
    return rcpp_result_gen;
END_RCPP
}
// trit2vec
IntegerVector trit2vec(int t);
RcppExport SEXP _MixedUpParents_trit2vec(SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(trit2vec(t));
    return rcpp_result_gen;
END_RCPP
}
// pgp_rcpp
List pgp_rcpp(IntegerMatrix IXG, NumericMatrix AF, IntegerVector isD, NumericMatrix AD, IntegerVector debug);
RcppExport SEXP _MixedUpParents_pgp_rcpp(SEXP IXGSEXP, SEXP AFSEXP, SEXP isDSEXP, SEXP ADSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type IXG(IXGSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type AF(AFSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type isD(isDSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type AD(ADSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(pgp_rcpp(IXG, AF, isD, AD, debug));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _MixedUpParents_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MixedUpParents_calc_HOT", (DL_FUNC) &_MixedUpParents_calc_HOT, 3},
    {"_MixedUpParents_intersect_ancestry_intervals_rcpp", (DL_FUNC) &_MixedUpParents_intersect_ancestry_intervals_rcpp, 9},
    {"_MixedUpParents_trit2vec", (DL_FUNC) &_MixedUpParents_trit2vec, 1},
    {"_MixedUpParents_pgp_rcpp", (DL_FUNC) &_MixedUpParents_pgp_rcpp, 5},
    {"_MixedUpParents_rcpp_hello", (DL_FUNC) &_MixedUpParents_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_MixedUpParents(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
