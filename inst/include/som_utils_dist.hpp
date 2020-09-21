#ifndef SOM_UTILS_DIST_HPP
#define SOM_UTILS_DIST_HPP


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <string>



namespace SOM_UTILS_DIST {

// *** Distances from one rowvec to another 
// Squared Euclidean 
inline double dist_L22(const arma::rowvec& x, const arma::rowvec& y, double bias = 0.0) {
  return arma::accu(arma::square(x - y)) - bias; 
}

// Euclidean
inline double dist_L2(const arma::rowvec& x, const arma::rowvec& y, double bias = 0.0) {
  double d2 = arma::accu(arma::square(x - y)); 
  return std::sqrt(d2) - bias; 
}

// Squared Mahalanobis
inline double dist_Mahalanobis2(const arma::rowvec& x, const arma::mat& S, const arma::rowvec& y, double bias = 0.0) {
  arma::rowvec z = x - y; 
  return arma::as_scalar( z * S * z.t() - bias );
}

// Mahalanobis 
inline double dist_Mahalanobis(const arma::rowvec& x, const arma::mat& S, const arma::rowvec& y, double bias = 0.0) {
  arma::rowvec z = x - y; 
  return std::sqrt(arma::as_scalar( z * S * z.t())) - bias;
}

// L1 - sum all elements
inline double dist_L1(const arma::rowvec& x, const arma::rowvec& y, double bias = 0.0) {
  return arma::accu( arma::abs(x-y) ) - bias;
}

// LInf - max all elements
inline double dist_LInf(const arma::rowvec& x, const arma::rowvec& y, double bias = 0.0) {
  return arma::max(arma::abs(x-y)) - bias;
}


// Parallel to find distances (make sure set with more rows is given as X1)
struct distmat_worker : public RcppParallel::Worker {
  
  const arma::mat& X1;
  const arma::mat& X2;
  std::string which_dist;
  
  // Possible inputs, if user requests to use  
  arma::vec bias;
  arma::rowvec X1min, X1max, X1rng, X2min, X2max, X2rng; 
  
  // Internal variables 
  unsigned int nX1, nX2, dim;
  bool use_bias, use_scaling; 
  
  // output container
  arma::mat DIST;
  
  
  // Constructor
  distmat_worker(const arma::mat& X1, const arma::mat& X2, std::string which_dist_) 
    : X1(X1), X2(X2), which_dist(which_dist_) {
    
    // Check that distance is valid 
    if(!(which_dist=="L2" || which_dist =="L22" || which_dist=="L1" || which_dist=="LInf")) Rcpp::stop("which_dist not valid.");
    
    // Check the X1 and X2 have same dimensions 
    dim = X1.n_cols; 
    if(X2.n_cols != dim) Rcpp::stop("ncol(X1) != ncol(X2)");
    
    // Initialize varibles
    nX1 = X1.n_rows; nX2 = X2.n_rows; 
    use_bias = false;
    use_scaling = false; 
    DIST.set_size(nX1,nX2);
  }
  
  // Set bias, if using 
  void set_bias(const arma::vec& bias_) {
    if(bias_.n_elem != nX2) Rcpp::stop("length(bias) != nrow(X2)");
    this->bias = bias_; 
    this->use_bias = true; 
  }
  
  // Set scaling ranges, if using
  void set_ranges(const arma::rowvec& X1min_, const arma::rowvec& X1max_, const arma::rowvec& X2min_, const arma::rowvec& X2max_) {
    
    // Set X1min 
    if(X1min_.n_elem == 1) {
      X1min.fill(X1min_[0]); 
    } else if(X1min_.n_elem == dim) {
      X1min = X1min_; 
    } else {
      Rcpp::stop("length(X1min) != ncol(X1)");
    }
    
    // Set X1max 
    if(X1max_.n_elem == 1) {
      X1max.fill(X1max_[0]); 
    } else if(X1max_.n_elem == dim) {
      X1max = X1max_; 
    } else {
      Rcpp::stop("length(X1max) != ncol(X1)");
    }
    
    // Set X1rng 
    X1rng = X1max - X1min; 
    
    // Set X2min 
    if(X2min_.n_elem == 1) {
      X2min.fill(X2min_[0]); 
    } else if(X2min_.n_elem == dim) {
      X2min = X2min_; 
    } else {
      Rcpp::stop("length(X2min) != ncol(X2)");
    }
    
    // Set X2max
    if(X2max_.n_elem == 1) {
      X2max.fill(X2max_[0]); 
    } else if(X2max_.n_elem == dim) {
      X2max = X2max_; 
    } else {
      Rcpp::stop("length(X2max) != ncol(X2)");
    }
    
    // Set Wrng 
    X2rng = X2max - X2min; 
    
    this->use_scaling = true; 
    
    return; 
  }
  
  // Scale a single data vector from X1 to the X2 range 
  arma::rowvec scale_x1_to_x2(unsigned int i) {
    arma::rowvec x = (X1.row(i) - X1min) / X1rng % X2rng + X2min; 
    return x; 
  }
  

  // Compute distance for a single (i,j) entry of the distance matrix 
  void dist_ij(unsigned int i, unsigned int j) {
    
    // Strip out the row vector X1(i). 
    // If we are scaling, do it 
    arma::rowvec x1i; 
    if(use_scaling) {
      x1i = scale_x1_to_x2(i);
    } else {
      x1i = X1.row(i); 
    }
    
    // Set the j-th bias to 0. If we are using bias, overwrite with the actual value. 
    double biasj = 0.0; 
    if(use_bias) biasj = bias(j);  
    
    // Now juse determine which distance we need to compute 
    if(which_dist == "L22") {
      DIST(i,j) = dist_L22(x1i, X2.row(j), biasj);
    } else if(which_dist=="L2") {
      DIST(i,j) = dist_L2(x1i, X2.row(j), biasj);
    } else if(which_dist == "L1") {
      DIST(i,j) = dist_L1(x1i, X2.row(j), biasj);
    } else if(which_dist == "LInf") {
      DIST(i,j) = dist_LInf(x1i, X2.row(j), biasj);
    }
    
  }
  

  // Parallel operator 
  void operator()(std::size_t begin, std::size_t end) {
    unsigned int i,j; 
    for(unsigned int k=begin; k<end; ++k) {
      i = k % nX1;
      j = k / nX1;
      
      dist_ij(i,j); 
    }
  }
  
  // Parallel invoker 
  void calc_parallel() {
    RcppParallel::parallelFor(0, nX1*nX2, *this);
  }
  
  // Serial invoker 
  void calc_serial() {
    unsigned int i,j; 
    for(unsigned int k=0; k<nX1*nX2; ++k) {
      i = k % nX1;
      j = k / nX1;
      
      dist_ij(i,j); 
    }
  }
  
};


// Geodesic graph distance
inline arma::umat geodesicdist(const arma::umat& ADJ, bool weighted = false, bool directed = false) {
  // Calculate the geodesic distances between vertices, given an adjacency matrix 
  
  // Tap into igraph package, setup functions 
  Rcpp::Environment igpkg = Rcpp::Environment::namespace_env("igraph");
  Rcpp::Function igfxn_graph_adjacency = igpkg["graph.adjacency"];
  Rcpp::Function igfxn_distances = igpkg["distances"];
  
  
  
  // Define the graph with given adjacency
  Rcpp::List ig; 
  if(weighted && directed) {
    ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","directed"), Rcpp::Named("diag", false), Rcpp::Named("weighted",true));  
  } else if(weighted && !directed) {
    ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","undirected"), Rcpp::Named("diag", false), Rcpp::Named("weighted",true));  
  } else if(!weighted && directed) {
    ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","directed"), Rcpp::Named("diag", false), Rcpp::Named("weighted",R_NilValue));  
  } else {
    ig = igfxn_graph_adjacency(Rcpp::Named("adjmatrix", ADJ), Rcpp::Named("mode","undirected"), Rcpp::Named("diag", false), Rcpp::Named("weighted",R_NilValue));  
  }
  
  
  // Compute shortest paths between all vertices 
  Rcpp::IntegerMatrix SPDISTRCPP = igfxn_distances(Rcpp::Named("graph", ig)); 
  
  // Convert to umat and return 
  arma::umat SPDIST = Rcpp::as<arma::umat>(SPDISTRCPP); 
  
  return SPDIST; 
}


} // close namespace

#endif