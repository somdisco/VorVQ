#include "vor_utils_gabriel.hpp"


//' Compute the Gabriel graph adjacency 
//' 
//' @param W the prototype matrix, prototype vectors in rows 
//' @param parallel boolean, whether to compute in parallel. Default = TRUE. 
//' 
//' @details 
//' The Gabriel graph is a proximity graph of points that is known to be a sub-graph 
//' of the Delaunay triangulation of the same set of points. 
//' 
//' @return a (binary, symmetric) adjacency matrix of dimension \code{nrow(W) x nrow(W)}.  Element 
//' \code{(i,j) = 1} IFF the prototype vectors \code{W[i,]} and \code{W[j,]} are Gabriel adjacent.  
//' Otherwise, entries are 0.  
//' 
//' @references
//' \insertRef{Gabriel1969}{VorVQ}
//' @export
// [[Rcpp::export]]
arma::umat Gabriel_ADJ(const arma::mat& W, bool parallel = true) {
  
  VOR_GABRIEL::Gabriel_ADJ_prlwkr wkr(W, 0); 
  if(parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  return wkr.GADJ; 
}