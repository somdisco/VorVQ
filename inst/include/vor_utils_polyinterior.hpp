#ifndef VORVQ_UTILS_POLYINTERIOR_HPP
#define VORVQ_UTILS_POLYINTERIOR_HPP

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]



#include "vor_utils_parallel.hpp"
#include "vor_utils_polydef.hpp"
#include "vor_utils_gurobi.hpp"




namespace VOR_POLYINTERIOR {

// ***** Check that points are interior 
inline bool is_interior_point(const arma::mat& A, const arma::vec& b, const arma::vec& x0) {
  // Compute slacks 
  arma::vec slacks = b - A * x0; 
  
  if(arma::min(slacks) > 0) return true;
  
  return false; 
}

// Parallel worker to check interior points for each First-Order Voronoi cell
struct vor1_check_interior_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& polydef_ADJ;
  const arma::uvec& i_list; // i_list denotes the idx (in W) of the Vor.cell we compute interior points for 
  const arma::vec& xlb; const arma::vec& xub;
  const arma::mat& X0; 
  
  // Internal variables
  unsigned int nprocess, d;

  // Outputs
  arma::uvec is_interior; 
  
  vor1_check_interior_prlwkr(const arma::mat& W_, const arma::umat& polydef_ADJ,
                             const arma::uvec& i_list_, const arma::vec& xlb_, const arma::vec& xub_, const arma::mat& X0) :
    W(W_), polydef_ADJ(polydef_ADJ), i_list(i_list_), xlb(xlb_), xub(xub_), X0(X0) {
    
    // Set sizes & perform additional checks
    nprocess = i_list.size(), d = W.n_cols;
    
    // Setup working container
    is_interior.set_size(nprocess); 
    is_interior.zeros(); 
    
  }
  
  
  // Function to compute the interior point for a single Vor1 Cell 
  // defined by i_list[i]
  void check_single_vor1_ip(unsigned int i) {
    
    // First obtain the polytope definition of this 1st-Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor1_polytope(W, polydef_ADJ, i_list(i), xlb, xub, A, b, cid);
    
    if(is_interior_point(A, b, X0.row(i).t())) is_interior[i] = 1; 
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = begin; i < end; i++) {
      
      check_single_vor1_ip(i);
      
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
    Rcpp::Rcout << "Checking interior points for vor1 cells ..." << std::endl; 
    
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
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Checking interior points for vor1 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      check_single_vor1_ip(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i_list(i)+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close solver function
};

// Parallel worker to check interior points for each Second-Order Voronoi cell
struct vor2_check_interior_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& polydef_ADJ;
  const arma::umat& ij_list; // rows of ij_list denote the idxs (in W) of the 2nd order Vor.cell we compute interior points for 
  const arma::vec& xlb; const arma::vec& xub;
  const arma::mat& X0; 
  
  // Internal variables
  unsigned int nprocess, d;
  
  // output 
  arma::uvec is_interior; 
  
  
  vor2_check_interior_prlwkr(const arma::mat& W_, const arma::umat& polydef_ADJ,
                             const arma::umat& ij_list_, const arma::vec& xlb_, const arma::vec& xub_, const arma::mat& X0) :
    W(W_), polydef_ADJ(polydef_ADJ), ij_list(ij_list_), xlb(xlb_), xub(xub_), X0(X0) {
    
    // Set sizes & perform additional checks
    nprocess = ij_list.n_rows, d = W.n_cols;
    if(nprocess != X0.n_rows) Rcpp::stop("nrow(ij_list) != nrow(X0)");
    
    is_interior.set_size(nprocess); 
    is_interior.zeros(); 
  }
  
  
  // Function to compute the interior point for a single Vor2 Cell 
  // defined by ij_list[i,]
  void check_single_vor2_ip(unsigned int i) {
    
    // First obtain the polytope definition of this 2nd Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor2_polytope(W, polydef_ADJ, ij_list(i,0), ij_list(i,1), xlb, xub, A, b, cid);
    
    if(is_interior_point(A, b, X0.row(i).t())) is_interior[i] = 1; 
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = begin; i < end; i++) {
      
      check_single_vor2_ip(i);
      
    } // close loop vor2 cells
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
    Rcpp::Rcout << "Checking interior points for vor2 cells ..." << std::endl; 
    
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
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Checking interior points for vor2 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      check_single_vor2_ip(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close solver function
  
  
};



// ***** Chebyshev Centers 
// www.math.uwaterloo.ca/~hwolkowi/henry/teaching/w10/367.w10/367miscfiles/Lecture1.pdf
inline arma::vec cpp_polytope_chebyshev_center(const GRBEnv* env, const arma::mat& A, const arma::vec& b, bool verbose) {
  
  // Strip out dimensions of polytope
  int M = A.n_rows, N = A.n_cols;
  // Initialize output 
  arma::vec center(N); center.zeros();
  
  
  try{
    
    // Initialize a Gurobi model
    GRBModel model = GRBModel(*env);
    if(!verbose) model.getEnv().set(GRB_IntParam_OutputFlag, 0); // Turn off
    std::vector<GRBVar> x(N+1);
    
    // Add the "y" variables, set all objective coefficients to 0 to start
    for(int n=0; n<N; ++n) {
      x[n] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x"+std::to_string(n)); //
    }
    x[N] = model.addVar(0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "x"+std::to_string(N)); //
    
    // Add the constraints 
    for(int m=0; m<M; ++m) {
      // Add the linear constraints
      GRBLinExpr constr = 0.0; // Initialize this linear constrant Aix <= bi   
      for(int n=0; n<N; ++n) {
        constr += A(m,n) * x[n]; // Fill up constraint i
      }
      
      double norm_Am = std::sqrt(arma::accu(arma::square(A.row(m)))); 
      constr += x[N] * norm_Am;
      
      // Add the constraint we built above
      model.addConstr(constr, GRB_LESS_EQUAL, b(m), "c"+std::to_string(m));
    }
    
    
    // The objective is to maximize each slack 
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    
    // Solve the model
    model.optimize();
    
    // Check if model wasn't solved.  If not, continue
    if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      for(int n=0; n<N; ++n) {
        center(n) += x[n].get(GRB_DoubleAttr_X);
      }
    }
    
  }
  catch (GRBException e)
  {
    if(verbose) {
      Rcpp::Rcout << "Error code = " << e.getErrorCode() << std::endl;
      Rcpp::Rcout << e.getMessage() << std::endl;
    }
  }
  catch (...)
  {
    if(verbose) {
      Rcpp::Rcout << "Exception during optimization" << std::endl;
    }
  }
  
  
  // Return the center.
  return center;    
}

// Parallel worker to compute Chebyshev Centers for each First-Order Voronoi cell
struct vor1_chebyshev_center_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& polydef_ADJ; // the adjacency matrix controlling definition of polytopes. Can be empty. 
  const arma::uvec& i_list; // i_list denotes the idx (in W) of the Vor.cell we compute interior points for 
  const arma::vec& xlb; const arma::vec& xub;
  
  // Internal variables
  unsigned int nprocess, d;
  //GRBEnv* env;
  
  // Outputs
  arma::mat C; // matrix of analytic centers, one row per element of i_list 
  
  vor1_chebyshev_center_prlwkr(const arma::mat& W_, const arma::umat& polydef_ADJ,
                             const arma::uvec& i_list_, const arma::vec& xlb_, const arma::vec& xub_) :
    W(W_), polydef_ADJ(polydef_ADJ), i_list(i_list_), xlb(xlb_), xub(xub_) {
    
    // Set sizes & perform additional checks
    nprocess = i_list.size(), d = W.n_cols;
    
    // Setup working container
    C.set_size(nprocess, d);
    
    // Initialize a Gurobi environment
    // FILE* mystdout = stdout;
    // stdout = fopen("/dev/null", "w");
    // std::fclose(stdout);
    // env = new GRBEnv();
    // stdout = mystdout;
  }
  
  
  // Function to compute the interior point for a single Vor1 Cell 
  // defined by i_list[i]
  void single_vor1_cheby(unsigned int i) {
    
    // First obtain the polytope definition of this 1st-Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor1_polytope(W, polydef_ADJ, i_list(i), xlb, xub, A, b, cid);  

    
    // Determine the interior point 
    arma::vec tmp = cpp_polytope_chebyshev_center(GUROBIENV, A, b, false);
    
    // Store it 
    C.row(i) = tmp.t(); 
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = begin; i < end; i++) {
      
      single_vor1_cheby(i);
    
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
    Rcpp::Rcout << "Calculating Chebyshev centers for vor1 cells ..." << std::endl; 
    
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
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating Chebyshev centers for vor1 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      single_vor1_cheby(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i_list(i)+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close solver function
};

// Parallel worker to compute Chebyshev Centers for each Second-Order Voronoi cell
struct vor2_chebyshev_center_prlwkr : RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W;
  const arma::umat& polydef_ADJ;
  const arma::umat& ij_list; // rows of ij_list denote the idxs (in W) of the 2nd order Vor.cell we compute interior points for 
  const arma::vec& xlb; const arma::vec& xub;
  
  // Internal variables
  unsigned int nprocess, d;
  //GRBEnv* env;
  
  
  // Outputs
  arma::mat C; // matrix of analytic centers, one row per element of i_list 
  
  vor2_chebyshev_center_prlwkr(const arma::mat& W_, const arma::umat& polydef_ADJ,
                             const arma::umat& ij_list_, const arma::vec& xlb_, const arma::vec& xub_) :
    W(W_), polydef_ADJ(polydef_ADJ), ij_list(ij_list_), xlb(xlb_), xub(xub_) {
    
    // Set sizes & perform additional checks
    nprocess = ij_list.n_rows, d = W.n_cols;
    
    // Setup working container
    C.set_size(nprocess, d);
    
    // Initialize a Gurobi environment
    // FILE* mystdout = stdout;
    // stdout = fopen("/dev/null", "w");
    // std::fclose(stdout);
    // env = new GRBEnv();
    // stdout = mystdout;
  }
  
  
  // Function to compute the interior point for a single Vor2 Cell 
  // defined by ij_list[i,]
  void single_vor2_cheby(unsigned int i) {
    
    // First obtain the polytope definition of this 2nd Order cell, with Delaunay restriction & global bounds
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor2_polytope(W, polydef_ADJ, ij_list(i,0), ij_list(i,1), xlb, xub, A, b, cid);
    
    // Determine the interior point 
    arma::vec tmp = cpp_polytope_chebyshev_center(GUROBIENV, A, b, false);
    
    // Store it 
    C.row(i) = tmp.t(); 
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = begin; i < end; i++) {
      
      single_vor2_cheby(i);
      
    } // close loop vor2 cells
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
    Rcpp::Rcout << "Calculating Chebyshev centers for vor2 cells ..." << std::endl; 
    
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
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  void calc_serial() {
    
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating Chebyshev centers for vor2 cells ..." << std::endl; 
    
    unsigned int chunk_counter = 0;
    for(int i=0; i<nprocess; ++i) {
      
      // Check for abort 
      Rcpp::checkUserInterrupt(); 
      
      single_vor2_cheby(i);
      chunk_counter++; 
      
      //Rcpp::Rcout << "success" << std::endl;
      Rcpp::Rcout << i+1; 
      
      if(i == (nprocess-1) || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    } // close loop over cells
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    
  } // close solver function
  
};


}

#endif



