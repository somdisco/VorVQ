#include "vor_utils_polyinterior.hpp"



//' Chebyshev center of a polytope
//'
//' @description 
//' For a feasible set of linear equations Ax<=b that define a bounded polytope P, 
//' the Chebyshev center is the center the largest inscribed ball contained in P.
//' Polytope definition should be given as Ax <= b. 
//'
//' @param A the LHS matrix of polytope constraints
//' @param b the RHS vector of polytope constraints 
//' @return a vector (of the same dimension as A) giving the interior point 
//' @export 
//[[Rcpp::export]]
arma::vec polytope_Chebyshev_center(const arma::mat& A, const arma::vec& b) {
  
  // Call worker function
  arma::vec out = VOR_POLYINTERIOR::cpp_polytope_chebyshev_center(GUROBIENV, A, b, false);
  
  return out;
}


//' Check if point is interior
//' 
//' @description
//' For a feasible set of linear equations Ax<=b that define a bounded polytope P 
//' and a point x0, x0 is interior IFF b - Ax > 0.  
//' 
//' @param A the LHS matrix of polytope constraints
//' @param b the RHS vector of polytope constraints 
//' @param x the query point 
//' @return boolean = T if the query point satisfies the above condition 
//' @export 
//[[Rcpp::export]]
bool is_interior_point(const arma::mat& A, const arma::vec& b, const arma::vec& x) {
  
  return VOR_POLYINTERIOR::is_interior_point(A, b, x); 
  
}



/*
 //' Analytic center of a polytope
 //'
 //' @description 
 //' For a feasible set of linear equations Ax<=b that define a bounded polytope, 
 //' the analytic center is the feasible point that maximizes prod(b-Ax).
 //' Polytope definition should be given as Ax <= b. 
 //'
 //' @param A the LHS matrix of polytope constraints
 //' @param b the RHS vector of polytope constraints 
 //' @param x0 a starting point for the optimization, must be interior to the polytope defined by A and b 
 //' @return a vector (of the same dimension as A) giving the interior point 
 //' @export 
 //[[Rcpp::export]]
 arma::vec polytope_analytic_center(const arma::mat& A, const arma::vec& b, arma::vec x0) {
 
 // Call worker function
 bool success = VOR_POLYINTERIOR::cpp_polytope_analytic_center(x0, A, b, "newton");
 
 if(success) return x0;   
 
 return x0;
 }
*/
 
