#include "vor_utils_mvie.hpp"


//' Max volume inscribed ellipsoid of polytope
//'
//' Returns the matrix E and vector c such that the ellipsoid \eqn{(y - c)'E^-1 (y - c) <= 1} is inscribed in the polytope defined by Ax <= b with maximum volume.
//' Code taken from Yin Zhang's \verb{MVE} routine.
//'
//' @param A matrix of LHS constraints in H-representation of polytope
//' @param b vector of RHS constraints in H-representation of polytope
//' @param c0 a point inside the polytope (not on the boundary)
//' @param maxiter integer number of iterations to perform, default = 100
//' @param tol tolerance for Zhang's MVE code, see paper, default = 1e-4
//' @param fix_c0 boolean, whether to fix the given center at c0 during iteration, default = F
//' @param verbose boolean, whether to display the solver status at each iteration, default = F
//' @return a list with components
//' \itemize{
//' \item{E}{a matrix defining the rotation and scale of the maximum volume ellipsoid}
//' \item{c}{a vector defining the center of the maximum volume ellipsoid}
//' \item{status}{a flag returned from Zhang's solver code.  See details.}
//' }
//' @details \verb{maxiter} and \verb{tol} are control (convergence) parameters
//' for Yin Zhang's \verb{MVE} routine.
//'
//' The polytope defined by A & b must be bounded.
//'
//' status = 0 indicates successful convergence; = -1 means maxiter was reached before convergence.
//'
//' @references 
//' \insertRef{ZhangGao2003}{VorVQ}
//' @export
// [[Rcpp::export]]
Rcpp::List max_vol_inscr_ell(arma::mat A, arma::vec b, const arma::vec& c0, int maxiter = 100, double tol = 1e-4, bool fix_c0 = false, bool verbose = false) {
  
  // Initialize output containers
  arma::mat E;
  arma::vec c;
  
  // Call internal function 
  int status = VOR_MVIE::cpp_max_vol_inscr_ell(E, c, A, b, c0, maxiter, tol, fix_c0, verbose);
  
  // Format outputs
  Rcpp::List out; 
  out["E"] = E; 
  out["c"] = c; 
  out["status"] = status; 
  
  return out; 
}
