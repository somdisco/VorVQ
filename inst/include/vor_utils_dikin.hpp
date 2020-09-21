#ifndef VORVQ_UTILS_DIKIN_HPP
#define VORVQ_UTILS_DIKIN_HPP


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include "vor_utils_polydef.hpp"
#include "vor_utils_parallel.hpp"



namespace VOR_DIKIN {

inline arma::mat cpp_dikin_ell(const arma::mat& A, const arma::vec& b, const arma::vec& x0) {
  
  // Initialize output 
  arma::mat out;
  
  // Compute slacks
  arma::vec slacks = b - A * x0; 
  
  // If any slacks are <= 0, remove them 
  arma::uvec keep_these = arma::find(slacks >= 1e-10);
  if(keep_these.size() != A.n_rows) {
    
    arma::mat Armv = A.rows(keep_these);
    arma::mat brmv = b.elem(keep_these);
    // Compute slacks
    slacks = slacks.elem(keep_these); 
    // Square them
    slacks = arma::pow(slacks, 2.0);
    // Invert them 
    slacks = 1.0 / slacks; 
    // Form Dikin ellipsoid
    out = Armv.t() * arma::diagmat(slacks) * Armv; 
    // Make sure it's symmetric 
    out += out.t();
    out *= 0.5; 
    return arma::inv(out);
  }
  

  // Otherwise, form the matrix as normal 
  
  // Square the slacks 
  slacks = arma::pow(slacks, 2.0);
  // Invert them 
  slacks = 1.0 / slacks; 
  
  // Form Dikin ellipsoid
  out = A.t() * arma::diagmat(slacks) * A; 
  // Make sure it's symmetric 
  out += out.t();
  out *= 0.5; 
  
  return arma::inv(out);
}

// ***** Parallel worker to solve Max. Vol. Inscribed Ellipsoid problem for each 1st-Order Voronoi cell
struct vor1_dikin_ell_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& DADJ;
  const arma::vec& xlb; const arma::vec& xub;
  const arma::uvec& i_list; // i_list denotes the idx (in W) of the Vor.cell bounds found in LB & UB, assumed 0-based 
  const arma::mat& C;
  
  // Internal variables
  unsigned int nprocess, d;
  
  // Outputs
  arma::cube E; // the MVIE transformations 
  
  vor1_dikin_ell_prlwkr(const arma::mat& W, const arma::umat& DADJ,
                        const arma::vec& xlb, const arma::vec& xub,
                        const arma::uvec& i_list, const arma::mat& C) :
    W(W), DADJ(DADJ), xlb(xlb), xub(xub), i_list(i_list), C(C) {
    
    // Set sizes & perform additional checks
    nprocess = i_list.size(), d = W.n_cols;
    
    // Setup working containers
    E.set_size(d, d, nprocess);
  }
  
  
  
  // Function to compute the MVIE for a single Vor1 Cell 
  // defined by i_list[i]
  void single_vor1_dikin_ell(unsigned int i) {
    
    // First obtain the polytope definition of this 1st-Order cell, with Delaunay restriction & with global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor1_polytope(W, DADJ, i_list(i), xlb, xub, A, b, cid, false);
    
    // Compute the Dikin ellipsoid
    E.slice(i) = cpp_dikin_ell(A, b, C.row(i).t());
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = int(begin); i < int(end); i++) {
      
      single_vor1_dikin_ell(i);
      
    } // close loop vor1 cells
  } // close worker
  
  void calc_parallel_progress() {
    
    // nprocess = i_list.size(), set in constructor 
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    // Bust the set of indices in i_list up into chunks so we can report
    
    unsigned int chunk_size = nprocess / 10; 
    chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
    chunk_size = std::min(chunk_size, nprocess); 
    unsigned int chunk_counter = 0; 
    unsigned int ncompleted = 0; 
    Rcpp::Rcout << "Calculating Dikin ellipsoids for vor1 cells ..." << std::endl; 
    
    while(ncompleted < nprocess) {
      Rcpp::checkUserInterrupt(); 
      
      unsigned int from = chunk_counter * chunk_size; 
      unsigned int to = from + chunk_size; 
      to = std::min(to, nprocess); 
      
      RcppParallel::parallelFor(from, to, *this); 
      
      ncompleted += to - from; 
      
      chunk_counter++; 
      
      Rcpp::Rcout << int(double(ncompleted)/double(nprocess)*100) << "%";
      if(ncompleted == nprocess || chunk_counter % 5 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
    }
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  void calc_serial() {
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating Dikin ellipsoids for vor1 cells ..." << std::endl;
    for(unsigned int j=0; j<nprocess; ++j) {
      Rcpp::Rcout << j+1; 
      
      single_vor1_dikin_ell(j);
      
      if(j == (nprocess-1) || (j+1) % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    }
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  }
};


// ***** Parallel worker to solve Max. Vol. Inscribed Ellipsoid problem for each 2nd-Order Voronoi cell
struct vor2_dikin_ell_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& DADJ;
  const arma::vec& xlb; const arma::vec& xub;
  const arma::umat& ij_list; // assumed 0-based
  const arma::mat& C;
  
  // Internal variables
  unsigned int d, nprocess;
  
  // Outputs
  arma::cube E;
  
  vor2_dikin_ell_prlwkr(const arma::mat& W, const arma::umat& DADJ, 
                        const arma::vec& xlb, const arma::vec& xub,
                        const arma::umat& ij_list, const arma::mat& C) :
    W(W), DADJ(DADJ), xlb(xlb), xub(xub), ij_list(ij_list), C(C) {
    
    // Set sizes & perform additional checks
    nprocess = ij_list.n_rows, d = W.n_cols;
    
    // Setup working containers
    d = W.n_cols; nprocess = ij_list.n_rows;
    E.set_size(d, d, nprocess);
  }
  
  
  // Function to compute the MVIE for a single Vor2 Cell 
  void single_vor2_dikin_ell(unsigned int i) {
    
    // First obtain the polytope definition of this 2nd Order cell, with Delaunay restriction & without global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    //cpp_vor2_polytope(W, DADJ, ij_list(i,0), ij_list(i,1), A, b, cid, false);
    VOR_POLYDEF::cpp_vor2_polytope(W, DADJ, ij_list(i,0), ij_list(i,1), xlb, xub, A, b, cid, false);
    
    // Compute the Dikin ellipsoid 
    E.slice(i) = cpp_dikin_ell(A, b, C.row(i).t());
  }
  
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = int(begin); i < int(end); i++) {
      
      single_vor2_dikin_ell(i);
      
    } // close loop vor2 cells
  } // close operator
  
  void calc_parallel_progress() {
    
    // nprocess = ij_list.n_rows, set in constructor 
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    // Bust the set of indices in i_list up into chunks so we can report
    
    unsigned int chunk_size = nprocess / 10; 
    chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
    chunk_size = std::min(chunk_size, nprocess); 
    unsigned int chunk_counter = 0; 
    unsigned int ncompleted = 0; 
    Rcpp::Rcout << "Calculating Dikin ellipsoids for vor2 cells ..." << std::endl; 
    
    while(ncompleted < nprocess) {
      Rcpp::checkUserInterrupt(); 
      
      unsigned int from = chunk_counter * chunk_size; 
      unsigned int to = from + chunk_size; 
      to = std::min(to, nprocess); 
      
      RcppParallel::parallelFor(from, to, *this); 
      
      ncompleted += to - from; 
      
      chunk_counter++; 
      
      Rcpp::Rcout << int(double(ncompleted)/double(nprocess)*100) << "%";
      if(ncompleted == nprocess || chunk_counter % 5 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
    }
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  void calc_serial() {
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating Dikin ellipsoids for vor2 cells ..." << std::endl;
    for(unsigned int j=0; j<nprocess; ++j) {
      Rcpp::Rcout << j+1; 
      
      single_vor2_dikin_ell(j);
      
      if(j == (nprocess-1) || (j+1) % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    }
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  }
};

}

#endif