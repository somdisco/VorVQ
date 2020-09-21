// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/VorVQ_types.hpp"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// is_Delaunay_nhb
int is_Delaunay_nhb(const arma::mat& W, int iidx, int jidx, const arma::vec& lb, const arma::vec& ub, bool verbose);
RcppExport SEXP _VorVQ_is_Delaunay_nhb(SEXP WSEXP, SEXP iidxSEXP, SEXP jidxSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type iidx(iidxSEXP);
    Rcpp::traits::input_parameter< int >::type jidx(jidxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(is_Delaunay_nhb(W, iidx, jidx, lb, ub, verbose));
    return rcpp_result_gen;
END_RCPP
}
// Delaunay_ADJ
arma::umat Delaunay_ADJ(const arma::mat& W, const arma::vec& lb, const arma::vec& ub, bool parallel);
RcppExport SEXP _VorVQ_Delaunay_ADJ(SEXP WSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(Delaunay_ADJ(W, lb, ub, parallel));
    return rcpp_result_gen;
END_RCPP
}
// Dikin_ell
arma::mat Dikin_ell(const arma::mat& A, const arma::vec& b, const arma::vec& c);
RcppExport SEXP _VorVQ_Dikin_ell(SEXP ASEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(Dikin_ell(A, b, c));
    return rcpp_result_gen;
END_RCPP
}
// Gabriel_ADJ
arma::umat Gabriel_ADJ(const arma::mat& W, bool parallel);
RcppExport SEXP _VorVQ_Gabriel_ADJ(SEXP WSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(Gabriel_ADJ(W, parallel));
    return rcpp_result_gen;
END_RCPP
}
// max_vol_inscr_ell
Rcpp::List max_vol_inscr_ell(arma::mat A, arma::vec b, const arma::vec& c0, int maxiter, double tol, bool fix_c0, bool verbose);
RcppExport SEXP _VorVQ_max_vol_inscr_ell(SEXP ASEXP, SEXP bSEXP, SEXP c0SEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP fix_c0SEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_c0(fix_c0SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(max_vol_inscr_ell(A, b, c0, maxiter, tol, fix_c0, verbose));
    return rcpp_result_gen;
END_RCPP
}
// vor1_polytope
Rcpp::List vor1_polytope(const arma::mat& W, int iidx, Rcpp::Nullable<Rcpp::IntegerMatrix> DADJ, Rcpp::Nullable<Rcpp::NumericVector> lb, Rcpp::Nullable<Rcpp::NumericVector> ub, bool rmv_redundancies);
RcppExport SEXP _VorVQ_vor1_polytope(SEXP WSEXP, SEXP iidxSEXP, SEXP DADJSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rmv_redundanciesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type iidx(iidxSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerMatrix> >::type DADJ(DADJSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type rmv_redundancies(rmv_redundanciesSEXP);
    rcpp_result_gen = Rcpp::wrap(vor1_polytope(W, iidx, DADJ, lb, ub, rmv_redundancies));
    return rcpp_result_gen;
END_RCPP
}
// vor2_polytope
Rcpp::List vor2_polytope(const arma::mat& W, int iidx, int jidx, Rcpp::Nullable<Rcpp::IntegerMatrix> DADJ, Rcpp::Nullable<Rcpp::NumericVector> lb, Rcpp::Nullable<Rcpp::NumericVector> ub, bool rmv_redundancies);
RcppExport SEXP _VorVQ_vor2_polytope(SEXP WSEXP, SEXP iidxSEXP, SEXP jidxSEXP, SEXP DADJSEXP, SEXP lbSEXP, SEXP ubSEXP, SEXP rmv_redundanciesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type iidx(iidxSEXP);
    Rcpp::traits::input_parameter< int >::type jidx(jidxSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerMatrix> >::type DADJ(DADJSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< bool >::type rmv_redundancies(rmv_redundanciesSEXP);
    rcpp_result_gen = Rcpp::wrap(vor2_polytope(W, iidx, jidx, DADJ, lb, ub, rmv_redundancies));
    return rcpp_result_gen;
END_RCPP
}
// polytope_Chebyshev_center
arma::vec polytope_Chebyshev_center(const arma::mat& A, const arma::vec& b);
RcppExport SEXP _VorVQ_polytope_Chebyshev_center(SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(polytope_Chebyshev_center(A, b));
    return rcpp_result_gen;
END_RCPP
}
// is_interior_point
bool is_interior_point(const arma::mat& A, const arma::vec& b, const arma::vec& x);
RcppExport SEXP _VorVQ_is_interior_point(SEXP ASEXP, SEXP bSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_interior_point(A, b, x));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_vor_module();

static const R_CallMethodDef CallEntries[] = {
    {"_VorVQ_is_Delaunay_nhb", (DL_FUNC) &_VorVQ_is_Delaunay_nhb, 6},
    {"_VorVQ_Delaunay_ADJ", (DL_FUNC) &_VorVQ_Delaunay_ADJ, 4},
    {"_VorVQ_Dikin_ell", (DL_FUNC) &_VorVQ_Dikin_ell, 3},
    {"_VorVQ_Gabriel_ADJ", (DL_FUNC) &_VorVQ_Gabriel_ADJ, 2},
    {"_VorVQ_max_vol_inscr_ell", (DL_FUNC) &_VorVQ_max_vol_inscr_ell, 7},
    {"_VorVQ_vor1_polytope", (DL_FUNC) &_VorVQ_vor1_polytope, 6},
    {"_VorVQ_vor2_polytope", (DL_FUNC) &_VorVQ_vor2_polytope, 7},
    {"_VorVQ_polytope_Chebyshev_center", (DL_FUNC) &_VorVQ_polytope_Chebyshev_center, 2},
    {"_VorVQ_is_interior_point", (DL_FUNC) &_VorVQ_is_interior_point, 3},
    {"_rcpp_module_boot_vor_module", (DL_FUNC) &_rcpp_module_boot_vor_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_VorVQ(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}