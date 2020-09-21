#include "vor_utils_polydef.hpp"



//' First Order Voronoi Cell Polytope
//' Returns the polytope definition Ax <= b of a single Voronoi cell
//' @param W matrix of prototypes in rows
//' @param iidx Index of Voronoi cell, 1-based.
//' @param DADJ (optional) - Delaunay adjacency. If used, the H-representation of the returned polytope will be tight, only including constraints from Delaunay neighbors.
//' @param lb (optional) - Global lower bound of data. If given, the H-representation of the returned polytope will include this as a lower bound.
//' @param ub (optional) - same as xlb, but for upper bounds. Both must be given in order for the global bounds to be appended to the polytope representation.
//' @param rmv_redundancies boolean, whether to check the system for redundancies, and remove any found 
//' @return a list with components A, b, cid
//' @export
// [[Rcpp::export]]
Rcpp::List vor1_polytope(const arma::mat& W, int iidx, Rcpp::Nullable<Rcpp::IntegerMatrix> DADJ = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericVector> lb = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> ub = R_NilValue, 
                         bool rmv_redundancies = false) {
  // Initialize containers to store results
  arma::mat A; arma::vec b; arma::uvec cid;
  
  // Parse the inputs to determine which function to call
  // First check whether a Delaunay adjacency matrix was given
  bool DADJ_restricted = false;
  arma::umat arma_DADJ;
  if(DADJ.isNotNull()) {
    DADJ_restricted = true;
    // We must copy the input adjacency matrix to a new object to use it
    // because it was given as optionally NULL
    Rcpp::IntegerMatrix tmp_DADJ(DADJ);
    arma_DADJ = Rcpp::as<arma::umat>(tmp_DADJ);
  }
  
  // Next check for global bounds. Both must exist, and be the right size, before they can be used
  bool use_global_bounds = false;
  arma::vec arma_xlb, arma_xub;
  if(lb.isNotNull() && ub.isNotNull()) {
    use_global_bounds = true;
    // Copy both over
    Rcpp::NumericVector tmp_xlb(lb); Rcpp::NumericVector tmp_xub(ub);
    if(tmp_xlb.size() != tmp_xub.size()) Rcpp::stop("lb & ub must be the same size, if given.");
    arma_xlb = Rcpp::as<arma::vec>(tmp_xlb);
    arma_xub = Rcpp::as<arma::vec>(tmp_xub);
  }
  
  // Now call the correct version of the function
  // Note that we pass iidx-1 to each function, since passing to internal C++ function which is 0-based indexing
  if(!DADJ_restricted && !use_global_bounds) {
    VOR_POLYDEF::cpp_vor1_polytope(W, iidx-1, A, b, cid, rmv_redundancies);
  } else if(!DADJ_restricted && use_global_bounds) {
    //Rcpp::Rcout << "Global bounds will be appended." << std::endl;
    VOR_POLYDEF::cpp_vor1_polytope(W, iidx-1, arma_xlb, arma_xub, A, b, cid, rmv_redundancies);
  } else if(DADJ_restricted && !use_global_bounds) {
    //Rcpp::Rcout << "Restricting to Delaunay neighbors." << std::endl;
    VOR_POLYDEF::cpp_vor1_polytope(W, arma_DADJ, iidx-1, A, b, cid, rmv_redundancies);
  } else if(DADJ_restricted && use_global_bounds) {
    //Rcpp::Rcout << "Global bounds will be appended, Restricting to Delaunay neighbors." << std::endl;
    VOR_POLYDEF::cpp_vor1_polytope(W, arma_DADJ, iidx-1, arma_xlb, arma_xub, A, b, cid, rmv_redundancies);
  }
  
  
  Rcpp::List out;
  out["A"] = A;
  out["b"] = b;
  out["cid"] = cid+1; // +1 since the indices are returned from 0-based internal C++ function
  return out;
}


//' Second Order Voronoi Cell Polytope
//' Returns the polytope definition Ax <= b of a second order Voronoi cell
//' @param W matrix of prototypes in rows
//' @param DADJ Delaunay adjacency. The H-representation of the returned polytope will be tight, only including constraints from Delaunay neighbors.
//' @param iidx Index of 1st Voronoi cell
//' @param jidx Index of 2nd Voronoi cell
//' @param lb - Global lower bound of data. If given, the H-representation of the returned polytope will include this as a lower bound.
//' @param ub - same as xlb, but for upper bounds. Both must be given in order for the global bounds to be appended to the polytope representation.
//' @param rmv_redundancies boolean, whether to check the system for redundancies, and remove any found 
//' @return a list with components A, b, cid
//' @export
// [[Rcpp::export]]
Rcpp::List vor2_polytope(const arma::mat& W, int iidx, int jidx,
                         Rcpp::Nullable<Rcpp::IntegerMatrix> DADJ = R_NilValue,  
                         Rcpp::Nullable<Rcpp::NumericVector> lb = R_NilValue, Rcpp::Nullable<Rcpp::NumericVector> ub = R_NilValue, 
                         bool rmv_redundancies = false) {
  
  // *** Parse the inputs to determine which function to call
  // First check whether a Delaunay adjacency matrix was given
  arma::umat arma_DADJ = arma::ones<arma::umat>(W.n_rows, W.n_rows); 
  arma_DADJ.diag().zeros(); 
  if(DADJ.isNotNull()) {
    // We must copy the input adjacency matrix to a new object to use it
    // because it was given as optionally NULL
    Rcpp::IntegerMatrix tmp_DADJ(DADJ);
    arma_DADJ = Rcpp::as<arma::umat>(tmp_DADJ);
  }
  
  // Next check for global bounds. Both must exist, and be the right size, before they can be used
  bool use_global_bounds = false;
  arma::vec arma_xlb, arma_xub;
  if(lb.isNotNull() && ub.isNotNull()) {
    use_global_bounds = true;
    // Copy both over
    Rcpp::NumericVector tmp_xlb(lb); Rcpp::NumericVector tmp_xub(ub);
    if(tmp_xlb.size() != tmp_xub.size()) Rcpp::stop("lb & ub must be the same size, if given.");
    arma_xlb = Rcpp::as<arma::vec>(tmp_xlb);
    arma_xub = Rcpp::as<arma::vec>(tmp_xub);
  }
  
  
  // ** Check that iidx & jidx are Delaunay neighbors. If not, can't proceed.
  if(!(arma_DADJ(iidx-1,jidx-1) > 0)) Rcpp::stop("Prototypes iidx & jidx are not Delaunay adjacent.");
  
  // Containers for polytope definition
  arma::mat A; arma::vec b; arma::uvec cid;
  if(use_global_bounds) {
    VOR_POLYDEF::cpp_vor2_polytope(W, arma_DADJ, iidx-1, jidx-1, arma_xlb, arma_xub, A, b, cid, rmv_redundancies);  
  } else {
    VOR_POLYDEF::cpp_vor2_polytope(W, arma_DADJ, iidx-1, jidx-1, A, b, cid, rmv_redundancies);  
  }
  
  
  
  Rcpp::List out;
  out["A"] = A;
  out["b"] = b;
  out["cid"] = cid+1; // +1 since the indices are returned from 0-based internal C++ function
  return out;
}


