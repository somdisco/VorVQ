#include "vor_utils_dikin.hpp"



//' Find Dikin ellipsoids inside a polytope
//' 
//' @description 
//' Returns a matrix E such that the ellipsoid \eqn{(y - c)'E^-1 (y - c) <= 1} is inscribed in the polytope defined by Ax <= b, 
//' with curvature influenced by the local geometry of the constraints at \eqn{c}
//' 
//' @param A LHS of polytope constraints
//' @param b RHS of polytope constraints
//' @param c interior point where Dikin ellipse is centered
//' @return the rotation matrix of the Dikin ellipsoid centered at \code{c}. 
//' 
//' 
//' @references 
//' \insertRef{Dikin1967}{VorVQ}
//' \insertRef{Boyd2004}{VorVQ}
//' @export 
//[[Rcpp::export]]
arma::mat Dikin_ell(const arma::mat& A, const arma::vec& b, const arma::vec& c) {
  arma::mat E = VOR_DIKIN::cpp_dikin_ell(A, b, c);
  return E;
}

