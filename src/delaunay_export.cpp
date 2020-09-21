#include "vor_utils_delaunay.hpp"



//' Determine if two prototypes are Delaunay neighbors
//'
//' @param W prototype matrix, with prototype vectors in rows
//' @param iidx the first index of the neighbor pair to test
//' @param jidx the second index of the neighbor pair to test
//' @param lb a vector of the global (data range) lower bounds, by dimension
//' @param ub a vector of the global (data range) upper bounds, by dimension
//' @param verbose boolean, whether to display the LP solver status
//' @return 1 = prototypes iidx and jidx are neighbors, 0 = they are not
//' @details iidj and jidx should be valid row indices of W indexed from 1 (not 0). \cr
//' The global lower and upper bounds (xlb and xub) are needed to ensure the polyope definition is bounded.
//' Otherwise, the LP we need to solve here could be unbounded, which produces incorrect results.
//' 
//' @references
//' \insertRef{Delaunay1934}{VorVQ}
//' 
//' \insertRef{Agrell1993}{VorVQ}
//' @export
// [[Rcpp::export]]
int is_Delaunay_nhb(const arma::mat& W, int iidx, int jidx, const arma::vec& lb, const arma::vec& ub, bool verbose=true) {

  // Compute Delaunay adjacency
  int result = VOR_DELAUNAY::is_Delaunay_nhb(GUROBIENV, W, iidx-1, jidx-1, lb, ub, verbose); // -1 since this interfaces with R

  return result;
}






//' Build the Delaunay Adjacency Matrix of a set of prototypes
//'
//' @param W prototype matrix, with prototype vectors in rows
//' @param lb a vector of the global (data range) lower bounds, by dimension
//' @param ub a vector of the global (data range) upper bounds, by dimension
//' @param parallel whether to process in parallel. Default = TRUE. 
//' @return A (binary, symmetric) adjacency matrix, nrows = ncols = nrow(W), giving
//' the Delaunay adjacencies between all prototypes in W. \cr
//' The global lower and upper bounds (xlb and xub) are needed to ensure the polyope definition is bounded.
//' Otherwise, the LP we need to solve here could be unbounded, which produces incorrect results.
//' 
//' @references
//' \insertRef{Delaunay1934}{VorVQ}
//' 
//' \insertRef{Agrell1993}{VorVQ}
//' @export
// [[Rcpp::export]]
arma::umat Delaunay_ADJ(const arma::mat& W, const arma::vec& lb, const arma::vec& ub, bool parallel = true) {

  // Initialize parallel structure
  VOR_DELAUNAY::Delaunay_ADJ_prlwkr wkr(GUROBIENV, W, lb, ub);
   
  if(parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  return wkr.DADJ; 
}










/*
// // To determine neighbor list of generators of 2nd-order masked Voronoi cells
// struct masked_Delaunay2_nhb_prlwkr : RcppParallel::Worker {
//   // Inputs
//   GRBEnv* env;
//   const arma::mat& Z;           // the 2nd order centroids (generators, or Voronoi sites)
//   const arma::umat& ij_list;    // The (i,j) idxs which define the vor2 cells whose centers are in Z
//   arma::umat CADJ;              // The CADJ adjacency of vor1 cells
//   const arma::vec& xlb, xub;    // Global bounds on X
//
//   // Calculate parameters
//   int nZ;
//
//   // Output
//   sivec delnhb;
//
//   // Constructor
//   masked_Delaunay2_nhb_prlwkr(GRBEnv* env_, const arma::mat& Z_, const arma::uvec& ij_list_, arma::umat CADJ_, const arma::vec& xlb_, const arma::vec& xub_) :
//     Z(Z_), ij_list(ij_list_), xlb(xlb_), xub(xub_) {
//     env = env_; // Set the Gurobi environment
//
//     // Make sure that ij_list does not reference any vor1 cells that DADJ does not represent
//     if(arma::any(ij_list.col(0) >= CADJ.n_rows)) Rcpp::stop("ij_list contains out-of-range indices.");
//     if(arma::any(ij_list.col(1) >= CADJ.n_rows)) Rcpp::stop("ij_list contains out-of-range indices.");
//
//     CADJ = CADJ_;
//     CADJ.diag().fill(1); // Allow each vor1 cell to be a neighbor of itself, this makes parallel looping easier
//
//     // Store # of 2nd-order points
//     nZ = Z.n_rows;
//
//     // Initialize output container
//     delnhb.resize(nZ);
//     for(int ii=0; ii<nZ; ++ii) {
//       delnhb[ii].resize(0);
//     }
//   };
//
//   // Parallel operator
//   void operator()(std::size_t begin, std::size_t end) {
//
//     // Loop over the points of Z in this block
//     for(int ii = begin; ii < end; ii++) {
//
//       // Get the vor1 cell to which point Z(ii) belongs
//       int vor1i = iZ(ii);
//
//       // Loop over all other points Z(jj) for which jj > ii
//       for(int jj=(ii+1); jj<nZ; ++jj) {
//
//         // Get the vor1 cell to which Z(jj) belongs
//         int vor1j = iZ(jj);
//
//         // If vor1i & vor1j are not first-order Delaunay neighbors, then
//         // ii & jj cannot be second-order neighbors. Skip.
//         if(DADJ(vor1i,vor1j) == 0) continue;
//
//         // Otherwise, they can be 2nd-order neighbors.  Compute the LP
//         int is_nhb = cpp_is_Delaunay_nhb(env, Z, ii, jj, xlb, xub, false);
//
//         // If the LP returned 1, ii & jj are neighbors. Add the adjacency to the list
//         if(is_nhb == 1) {
//           delnhb[ii].push_back(jj);
//         }
//
//       } // close loop over neighbors jj
//
//     } // close loop over ii
//   } // close worker
//
//   // To "symmetrize" the neighbor list
//   // Only perform this task after parallel computation,
//   // which only returns the adjacencies in the upper triangle of the DADJ2 adj. matrix
//   void symmetrize_delnhb() {
//
//     for(int ii=0; ii<nZ; ++ii) {
//
//       for(int jj=0; jj<delnhb[ii].size(); ++jj) {
//         // The connection ii -> jj already exists in the list.
//         // Must add jj -> ii
//         delnhb[jj].push_back(ii);
//       } // close jj loop
//     } // close ii loop
//
//   }
//
//   // To sort the delnhb list
//   void sort_delnhb() {
//     for(int ii=0; ii<nZ; ++ii) {
//       std::sort(delnhb[ii].begin(), delnhb[ii].end());
//     }
//   }
//
// };
//
// // [[Rcpp::export]]
// sivec Delaunay2_ADJ(const arma::mat& Z, const arma::uvec& iZ, const arma::umat& DADJ, const arma::vec& xlb, const arma::vec& xub) {
//
//   // Since iZ is coming from 1-based R, convert it to 0-based to pass to c++
//   arma::uvec iZ_cpp = iZ;
//   iZ_cpp -= 1;
//
//   // Initialize a Gurobi model
//   GRBEnv* env = 0;
//   env = new GRBEnv();
//
//   // Initialize parallel struct
//   Delaunay2_nhb_prlwkr wkr(env, Z, iZ_cpp, DADJ, xlb, xub);
//
//   // Test for 2nd-order Delaunay adjacencies
//   RcppParallel::parallelFor(0, Z.n_rows, wkr);
//
//   // Symmetrize the adjacency list just computed
//   wkr.symmetrize_delnhb();
//
//   // Sort the adjacency list just computed
//   wkr.sort_delnhb();
//
//   // Extract the adjacency list just computed & convert its indices back to 1-based to return to R
//   sivec delnhb = wkr.delnhb;
//
//   for(int ii=0; ii<wkr.nZ; ++ii) {
//     for(int jj=0; jj<delnhb[ii].size(); ++jj) {
//       delnhb[ii][jj] += 1;
//     }
//   }
//
//   // Clean up Gurobi environment
//   delete env;
//
//   // Return Delaunay adjacency matrix
//   return delnhb;
// }
//
//
*/


