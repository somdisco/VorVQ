#ifndef VORVQ_UTILS_GABRIEL_HPP
#define VORVQ_UTILS_GABRIEL_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include "vor_utils_parallel.hpp"


namespace VOR_GABRIEL {

// Code for the k-Gabriel Graph taken from here:
// ftp://ftp.cs.brown.edu/pub/techreports/96/cs96-28.pdf
inline bool is_Gabriel_nhb(const arma::mat& W, unsigned int iidx, unsigned int jidx, unsigned int korder=0) {
  // W is the matrix of all Voronoi centers
  
  unsigned int nW = W.n_rows;
  
  // Compute the distance from Wi to Wj
  // and midpoint between the two
  double d_ij = std::sqrt(arma::accu(arma::square(W.row(iidx)-W.row(jidx))));
  arma::rowvec mid = (W.row(iidx) + W.row(jidx)) / 2.0;
  
  // Loop over all centers (except i and j)
  // and compute distance to their midpoint.
  // If any point is found inside the sphere, a Gabriel edge cannot exist (so break and return FALSE)
  unsigned int interior_count = 0; 
  for(unsigned int k=0; k<nW; ++k) {
    if(k == iidx || k == jidx) continue;
    
    // Compute dist from Wk to mid
    double d_kmid = std::sqrt(arma::accu(arma::square(W.row(k)-mid)));
    
    if(d_kmid <= d_ij/2.0) {
      interior_count++; 
    } 
    
    if(interior_count > korder) return false;
  }
  
  return true;
}

// ***** Compute the full Delaunay adjacency *****
struct Gabriel_ADJ_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  unsigned int korder; 
  
  // Calculate parameters
  unsigned int nW;
  
  // Output
  arma::umat GADJ;
  
  // Constructor
  Gabriel_ADJ_prlwkr(const arma::mat& W_, unsigned int korder_ = 0) : W(W_) {
    nW = W.n_rows;
    GADJ = arma::zeros<arma::umat>(nW, nW);
    korder = korder_; 
  };
  
  // Test for Delaunay neighbors of j
  void neighbors_of_j(unsigned int j) {
    
    for(unsigned int k=j+1; k<nW; ++k) {
      
      bool is_nhb = is_Gabriel_nhb(W, j, k, korder);
      
      if(is_nhb) {
        GADJ(j,k) = 1;
        GADJ(k,j) = 1;
      }
    }
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    // Loop over neighbors of prototype iidx
    for(int j = int(begin); j < int(end); j++) {
      neighbors_of_j(j);
    } // close loop over neighbors
  } // close worker
  
  void calc_serial() {
    // log start time
    auto starttime = std::chrono::high_resolution_clock::now();
    
    Rcpp::Rcout << "Calculating Gabriel graph edges ..." << std::endl;
    for(unsigned int j=0; j<nW; ++j) {
      Rcpp::Rcout << j+1;
      
      neighbors_of_j(j);
      
      if(j == (nW-1) || (j+1) % 8 == 0) {
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
  
  
  void calc_parallel_progress() {
    
    // log start time
    auto starttime = std::chrono::high_resolution_clock::now();
    
    // Bust the set of indices in i_list up into chunks so we can report
    unsigned int nprocess = nW;
    
    unsigned int chunk_size = nprocess / 10;
    chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
    chunk_size = std::min(chunk_size, nprocess);
    unsigned int chunk_counter = 0;
    unsigned int ncompleted = 0;
    Rcpp::Rcout << "Calculating Gabriel graph edges ..." << std::endl;
    
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
  
  
};

}

#endif

