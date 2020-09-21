#ifndef VORVQ_UTILS_MVIE_HPP
#define VORVQ_UTILS_MVIE_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include "vor_utils_polydef.hpp"
#include "vor_utils_parallel.hpp"

namespace VOR_MVIE {

// Internal function: Componentwise min/max functions, needed for the code below 
inline void pmax(double testval, arma::vec& v) {
  for(unsigned int i=0; i<v.size(); ++i) {
    if(v[i] < testval) v[i] = testval;
  }
}


// ***** MVIE - Translation of Zhang's code
// Outputs are E = rotation matrix of ellipsoid)
// c = centroid of ellipsoid 
inline int cpp_max_vol_inscr_ell(arma::mat& E, arma::vec& c, arma::mat A, arma::vec b, const arma::vec& c0, int maxiter = 100, double tol = 1e-4, bool fix_c0 = false, bool verbose = false) {
  
  /*** Strip out dimensions & set some constants ***/
  int m = A.n_rows; 
  int n = A.n_cols;
  double bnrm = arma::norm(b,2); 
  double minmu = 1.e-8, tau0 = .75;  // these values were fixed in Zhang's code, left them that way here
  
  
  
  /*** Variable declarations ***/
  
  // The main variables which solve system (3.9) from Zhang's paper 
  c.resize(n); 
  //arma::vec c(n); // I have changed his variables x & x0 to c & c0 for "center"
  arma::vec y(m); 
  arma::vec z(m);
  
  // Storage for the ellipsoidal transform E s.t. (x-c)' E^-1 (x-c) <= 1
  // Zhang calls this matrix E2 in his code 
  //arma::mat E(n,n);
  E.resize(n,n);
  
  // Initialize helper variables needed during solution 
  arma::mat Q(m,m);
  arma::vec h(m);
  arma::vec yz(m), yh(m);
  double gap, rmu; 
  arma::vec R1(n), R2(m), R3(m);
  double r1, r2, r3, res, det_val, det_sign;
  arma::mat hprime(m,m);
  arma::mat N(m,m);
  arma::mat AtNM2(n,m);
  arma::vec R3Dy(m), R23(m); 
  arma::vec dc(n), dy(m), dz(m), Adc(m);
  
  
  
  // *** Setup for looping *** 
  
  // Center the linear system at 0
  arma::vec bmAc = b - A * c0;
  if(arma::any(bmAc <= 0)) Rcpp::stop("c0 not interior.");
  A.each_col() /= bmAc; 
  b.ones();
  bmAc = b;
  
  // Reset the main variables 
  c.zeros(); 
  y.ones(); 
  z.zeros();
  
  // Derive the starting values for the main helper variables 
  E = arma::inv(A.t() * arma::diagmat(y) * A);
  Q = A * E * A.t();
  h = arma::sqrt(Q.diag());
  
  // Print header if verbose
  if(verbose) {
    Rprintf("\n  Residuals:   Primal     Dual    Duality  logdet(E)\n");
    Rprintf("  --------------------------------------------------\n");
  }
  
  // The variable status determines the exit conditions
  //  0 = successful convergence 
  // -1 = maxiter was reached. Leave at this value unless tolerance test is successful. 
  int status = -1; 
  
  
  
  /*** Begin Looping ***/
  for(int iter=1; iter<=maxiter; ++iter) {
    
    // Only do on 1st iteration 
    if(iter == 1) {
      double t = arma::min(bmAc / h); 
      y /= std::pow(t,2.0); 
      h *= t;
      z = bmAc-h;
      pmax(1.e-1, z);
      Q *= std::pow(t,2.0); 
    }
    
    
    yz = y % z; 
    yh = y % h;
    gap = arma::accu(yz) / double(m);
    rmu = std::min(0.5, gap) * gap;
    rmu = std::max(rmu, minmu);
    
    // Set R1,R2,R3 & their norms
    // Here, Zhang used capital R for the vector and lowercase r for its norm
    // In the paper, lowercase r denotes one of the auxiliary vectors in system 3.22
    if(fix_c0) {
      r1 = 0.0; // If not updating the center c (x), do not let it affect convergence criteria 
    } else {
      R1 = -A.t() * yh;  
      r1 = arma::norm(R1, "inf");
    }
    
    R2 = bmAc - h - z;
    r2 = norm(R2, "inf");
    
    R3 = rmu - yz;
    r3 = norm(R3, "inf");
    
    res = std::max(r1, std::max(r2, r3));
    arma::log_det(det_val, det_sign, E); 
    
    // Print iter status if verbose 
    if(verbose) {
      Rprintf("  iter %3i  ", iter);
      Rprintf("%9.1e %9.1e %9.1e  %9.3e\n", r2,r1,r3,det_val);
    }
    
    // Check for convergence 
    if(res < tol*(1 + bnrm) && rmu <= minmu) {
      if(verbose) Rprintf("  Converged!\n");
      status = 0; // 0 means successful convergence 
      break;
    }
    
    // *** Compute step directions *** 
    // h'(y), equation 3.14 from paper 
    hprime = Q % Q;
    hprime.each_col() /= h; 
    hprime *= -0.5; 
    
    // N(y), equation 3.17 from paper 
    N = hprime; 
    N.each_col() %= y;
    N.diag() += h; 
    
    // Compute inv(M_2), equation 3.23 from paper
    // Instead of making another m x m array to store in memory, re-use the container we made for hprime, 
    // as it is not needed for further calculation. 
    // From hereon in the code, hprime = inv(M_2)
    hprime *= -1.0; 
    hprime.diag() += z/y; 
    hprime = arma::inv(hprime);
    
    // Pre-compute the product A^T * N * inv(M_2)
    AtNM2 = A.t() * N * hprime; 
    
    // These quantities are needed several times in system 3.22 from paper, pre-compute them here 
    R3Dy = R3 / y; 
    R23 = R2 - R3Dy;
    
    // dx step (3.22a): If not updating center c=x, do not solve the system to save time 
    if(fix_c0) {
      dc.zeros(); 
      Adc.zeros();
    } else {
      dc = arma::solve(AtNM2 * A, R1 + AtNM2 * R23);
      Adc = A * dc; 
    }
    
    // dy & dz steps, 3.22b-c 
    dy = hprime * (Adc - R23);
    dz = R3Dy - z/y % dy; 
    
    
    // Compute step sizes
    double ac = -1.0 / std::min( arma::min(-Adc / bmAc), -0.5);
    double ay = -1.0 / std::min( arma::min(dy / y), -0.5);
    double az = -1.0 / std::min( arma::min(dz / z), -0.5);
    double tau = std::max(tau0, 1.0-res);
    double astep = tau * std::min(std::min(std::min(1.0, ac), ay), az);
    
    // Update 
    c += astep * dc;
    y += astep * dy;
    z += astep * dz;
    
    // Recompute slacks 
    bmAc = bmAc - astep*Adc; 
    
    // Update E, Q, h. 
    // E is really only computed to monitor its determinant, y is the real workhorse 
    E = arma::inv(A.t() * arma::diagmat(y) * A); // Want to use inv_sympd here, but it is too sensitive to slight departures from pos.semi.def.
    Q = A * E * A.t();
    h = arma::sqrt(Q.diag());
    
  }
  
  // We centered at 0 at the start, put back here 
  c += c0; 
  
  return status; 
}


// ***** Parallel worker to solve Max. Vol. Inscribed Ellipsoid problem for each 1st-Order Voronoi cell
struct vor1_mvie_solve_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& DADJ;
  const arma::uvec& i_list; // i_list denotes the idx (in W) of the Vor.cell bounds found in LB & UB
  const arma::vec& xlb; const arma::vec& xub;
  const arma::mat& C0;
  int maxiter;
  double tol;
  bool fix_c0;
  bool use_DADJ, use_bounds;
  
  // Internal variables
  unsigned int nprocess, d;
  
  // Outputs
  arma::cube E; // the MVIE transformations 
  arma::mat C; // the MVIE centers
  arma::ivec status; // solver status from MVIE code 
  arma::vec logdet; 
  arma::vec volratio; 
  
  
  vor1_mvie_solve_prlwkr(const arma::mat& W_, const arma::umat& DADJ_,
                         const arma::uvec& i_list_, const arma::vec& xlb_, const arma::vec& xub_,
                         const arma::mat& C0_,
                         int maxiter_ = 100, double tol_ = 1e-4, 
                         bool fix_c0_ = false, bool use_DADJ_ = true, bool use_bounds_ = true) :
    W(W_), DADJ(DADJ_), i_list(i_list_), xlb(xlb_), xub(xub_), C0(C0_) {
    
    use_DADJ = use_DADJ_;
    use_bounds = use_bounds_;
    
    // Set sizes & perform additional checks
    nprocess = i_list.size(), d = W.n_cols;
    
    // Set the mve parameters
    maxiter = maxiter_;
    tol = tol_;
    fix_c0 = fix_c0_;
    
    // Setup working containers
    status.set_size(nprocess); 
    C.set_size(nprocess, d);
    E.set_size(d, d, nprocess);
    logdet.set_size(nprocess); logdet.zeros(); 
    volratio.set_size(nprocess); volratio.zeros(); 
  }
  
  
  
  // Function to compute the MVIE for a single Vor1 Cell 
  // defined by i_list[i]
  void single_vor1_mvie(int i) {
    
    // First obtain the polytope definition of this 1st-Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    if(use_DADJ && use_bounds) {
      VOR_POLYDEF::cpp_vor1_polytope(W, DADJ, i_list(i), xlb, xub, A, b, cid, false);  
    } else if(!use_DADJ && use_bounds) {
      VOR_POLYDEF::cpp_vor1_polytope(W, i_list(i), xlb, xub, A, b, cid, false);
    } else if(!use_DADJ && !use_bounds) {
      VOR_POLYDEF::cpp_vor1_polytope(W, i_list(i), A, b, cid, false);
    }
    
    
    // Run the MVIE code & store the center & rotation matrix
    arma::mat tmpE; 
    arma::vec tmpc; 
    status[i] = cpp_max_vol_inscr_ell(tmpE, tmpc, A, b,  C0.row(i).t(), maxiter, tol, fix_c0, false);
    E.slice(i) = tmpE; 
    C.row(i) = tmpc.t();
    
    // Compute the log determinant of the rotation matrix 
    arma::vec singularvals; 
    bool svd_pass = arma::svd(singularvals, tmpE);
    if(svd_pass) logdet[i] = arma::accu(arma::log(singularvals)); else logdet[i] = arma::datum::nan; 
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(int i = int(begin); i < int(end); i++) {
      
      single_vor1_mvie(i);
      
    } // close loop vor1 cells
  } // close worker
  
  void calc_parallel_progress() {
    
    // nprocess = length(i_list), set during constructor
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    // initialize chunk variables 
    unsigned int chunk_size = nprocess / 10; 
    chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
    chunk_size = std::min(chunk_size, nprocess); 
    unsigned int chunk_counter = 0; 
    unsigned int ncompleted = 0; 
    Rcpp::Rcout << "Calculating max vol incribed ellipsoids for vor1 cells ..." << std::endl; 
    
    // loop over chunks
    while(ncompleted < nprocess) {
      // Check for abort 
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
    } // close while 
    
    // Compute volume ratio 
    this->set_volratio(); 
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating max vol incribed ellipsoids for vor1 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      single_vor1_mvie(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i_list(i)+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // Compute volume ratio 
    this->set_volratio(); 
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close serial solver
  
  void set_volratio() {
    arma::vec logvol = 0.5 * this->logdet; 
    this->volratio = arma::exp(logvol - logvol.max()); 
    this->volratio = this->volratio / arma::accu(this->volratio); 
  }
};

// ***** Parallel worker to solve Max. Vol. Inscribed Ellipsoid problem for each 2nd-Order Voronoi cell
struct vor2_mvie_solve_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& DADJ;
  const arma::umat& ij_list;
  const arma::vec& xlb; const arma::vec& xub;
  const arma::mat& C0;
  int maxiter;
  double tol;
  bool fix_c0;
  
  // Internal variables
  unsigned int d, nprocess;
  
  // Outputs
  arma::ivec status;
  arma::cube E;
  arma::mat C;
  arma::vec logdet; 
  arma::vec volratio; 
  
  vor2_mvie_solve_prlwkr(const arma::mat& W_, const arma::umat& DADJ_,
                         const arma::umat& ij_list_, const arma::vec& xlb_, const arma::vec& xub_,
                         const arma::mat& C0_,
                         int maxiter_ = 100, double tol_ = 1e-4, 
                         bool fix_c0_ = false) :
    W(W_), DADJ(DADJ_), ij_list(ij_list_), xlb(xlb_), xub(xub_), C0(C0_) {
    
    // Set sizes & perform additional checks
    nprocess = ij_list.n_rows, d = W.n_cols;
    
    // Set the mve parameters
    maxiter = maxiter_;
    tol = tol_; 
    
    fix_c0 = fix_c0_;
    
    // Setup working containers
    d = W.n_cols; nprocess = ij_list.n_rows;
    status.set_size(nprocess); status.zeros();
    C.set_size(nprocess, d);
    E.set_size(d, d, nprocess);
    logdet.set_size(nprocess); logdet.zeros(); 
    volratio.set_size(nprocess); volratio.zeros(); 
  }
  
  
  // Function to compute the MVIE for a single Vor2 Cell 
  void single_vor2_mvie(unsigned int i) {
    
    // First obtain the polytope definition of this 2nd Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor2_polytope(W, DADJ, ij_list(i,0), ij_list(i,1), xlb, xub, A, b, cid, false);
    
    // Run the MVIE code 
    // Use prototype W(i_list(i)) as a starting center point, since we know it's in the cell
    arma::mat tmpE; 
    arma::vec tmpc; 
    status[i] = cpp_max_vol_inscr_ell(tmpE, tmpc, A, b, C0.row(i).t(), maxiter, tol, fix_c0, false);
    
    E.slice(i) = tmpE; 
    C.row(i) = tmpc.t();
    
    arma::vec singularvals; 
    bool svd_pass = arma::svd(singularvals, tmpE);
    if(svd_pass) logdet[i] = arma::accu(arma::log(singularvals)); else logdet[i] = arma::datum::nan; 
  }
  
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(int i = int(begin); i < int(end); i++) {
      
      single_vor2_mvie(i);
      
    } // close loop vor2 cells
  } // close operator
  

  void calc_parallel_progress() {
    
    // nprocess = nrow(ij_list), set during constructor
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    // initialize chunk variables 
    unsigned int chunk_size; 
    if(nprocess < 1000) {
      chunk_size = nprocess / 10; 
    } else {
      chunk_size = nprocess / 20; 
    }
    chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
    chunk_size = std::min(chunk_size, nprocess); 
    unsigned int chunk_counter = 0; 
    unsigned int ncompleted = 0; 
    Rcpp::Rcout << "Calculating max vol incribed ellipsoids for vor2 cells ..." << std::endl; 
    
    // loop over chunks
    while(ncompleted < nprocess) {
      // Check for abort 
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
    } // close while 
    
    // Compute volume ratio 
    this->set_volratio();
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating max vol incribed ellipsoids for vor2 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      single_vor2_mvie(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // Compute volume ratio 
    this->set_volratio();
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close solver function
  
  void set_volratio() {
    arma::vec logvol = 0.5 * this->logdet; 
    this->volratio = arma::exp(logvol - logvol.max()); 
    this->volratio = this->volratio / arma::accu(this->volratio); 
  }
  
};




} // close namespace 


#endif

