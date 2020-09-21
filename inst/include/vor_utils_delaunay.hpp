#ifndef VORVQ_UTILS_DELAUNAY_HPP
#define VORVQ_UTILS_DELAUNAY_HPP


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]


#include "vor_utils_polydef.hpp"
#include "vor_utils_parallel.hpp"
#include "som_utils_dist.hpp"
#include "vor_utils_gurobi.hpp"


namespace VOR_DELAUNAY {

  //typedef std::vector< std::vector<int> > sivec;
  const arma::uvec emptyuvec = arma::uvec();
  
  
  // ***** Determine if two prototypes are Delaunay neighbors *****
  inline int is_Delaunay_nhb(const GRBEnv* env, const arma::mat& W, int iidx, int jidx, const arma::vec& lb, const arma::vec& ub, bool verbose) {
    
    // Get the constraints & RHS of the linear program
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor1_polytope(W, iidx, A, b, cid, false);
    //VOR_POLYDEF::vor1_polytope(A, b, cid, iidx, W, false); 
    int M = A.n_rows, N = A.n_cols;
    
    // Build the objective function coefficients for this pair
    arma::uvec coef_is_j = arma::find(cid == jidx);
    arma::vec coef = A.row(coef_is_j(0)).t();
    
    try{
      
      // Initialize a Gurobi model
      GRBModel model = GRBModel(*env);
      if(!verbose) model.getEnv().set(GRB_IntParam_OutputFlag, 0); // Turn off
      std::vector<GRBVar> x(N);
      
      // Add the variables
      for(int n=0; n<N; ++n) {
        x[n] = model.addVar(lb(n), ub(n), coef(n), GRB_CONTINUOUS, "x"+std::to_string(n));
      }
      
      // The objective is to minimize the costs
      model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
      
      // Add the constraints
      for(int m=0; m<M; ++m) {
        GRBLinExpr ntot = 0;
        for(int n=0; n<N; ++n) {
          ntot += A(m,n) * x[n];
        }
        model.addConstr(ntot, GRB_LESS_EQUAL, b(m), "c"+std::to_string(cid(m)));
      }
      
      // Optimize
      model.optimize();
      
      // Check if model wasn't solved.  If not, return 0 & exit
      if(model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        return 0;
      }
      
      // Extract the constraint
      GRBConstr constr = model.getConstrByName("c"+std::to_string(jidx));
      double constr_val = constr.get(GRB_DoubleAttr_Pi);
      
      // Return with 1 if this is a Delaunay neighbor
      if(constr_val > 0) {
        return 1;
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
    
    return 0;
    
  };
  
  
  // ***** Determine all Delaunay neighbors of prototype i *****
  // This only tests neighbors j > i. Must symmetrize the output to get actual DADJ 
  inline arma::uvec Delaunay_nhbs_of_i(const GRBEnv* env, const arma::mat& W, unsigned int iidx, const arma::vec& lb, const arma::vec& ub, bool verbose, arma::uvec jidxs = emptyuvec) {
    
    // If iidx = the last prototype, then there are no neighbors to test. just return an empty vector 
    if(iidx == (W.n_rows-1)) return arma::zeros<arma::uvec>(0); 
    
    // The input jidxs (possibly) contains a list of prototype indiced to restrict 
    // the Delauany testing to. e.g.: if given, only the adjacencies between 
    // (iidx, jidxs[0]), (iidx, jidxs[1]), etc will be tested. 
    // If not given, just populate the vector with all of the indices > iidx as default. 
    if(jidxs.n_elem == 0) {
      jidxs = arma::regspace<arma::uvec>(iidx+1, W.n_rows-1); 
    }
    
    // Get the constraints & RHS of the linear program
    arma::mat A; arma::vec b; arma::uvec cid;
    VOR_POLYDEF::cpp_vor1_polytope(W, iidx, A, b, cid, false);
    int M = A.n_rows, N = A.n_cols;
    
    // Empty output vector to start 
    std::vector<unsigned int> nhblist; 
    
    try{
      // Initialize a Gurobi model
      GRBModel model = GRBModel(*env);
      if(!verbose) model.getEnv().set(GRB_IntParam_OutputFlag, 0); // Turn off
      // The objective is to minimize the costs
      model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
      
      
      // Add the variables, initially with 1.0 objective coefficients
      std::vector<GRBVar> x(N);
      for(int n=0; n<N; ++n) {
        x[n] = model.addVar(lb(n), ub(n), 1.0, GRB_CONTINUOUS, "x"+std::to_string(n));
      }
      
      // Add the constraints
      for(int m=0; m<M; ++m) {
        GRBLinExpr ntot = 0;
        for(int n=0; n<N; ++n) {
          ntot += A(m,n) * x[n];
        }
        model.addConstr(ntot, GRB_LESS_EQUAL, b(m), "c"+std::to_string(cid(m)));
      }
      
      model.update(); 
      
      // Loop over all possible neighbors 
      for(unsigned int j=0; j<jidxs.n_elem; ++j) {
      //for(unsigned int j=(iidx+1); j<W.n_rows; ++j) {
        
        // Build the objective function coefficients for this pair
        //arma::uvec coef_is_j = arma::find(cid == j);
        arma::uvec coef_is_j = arma::find(cid == jidxs[j]);
        
        // Change the objective coefficients of the variables
        for(int n=0; n<N; ++n) {
          x[n].set(GRB_DoubleAttr_Obj, A(coef_is_j(0),n)); 
        }
        
        // Update and solve model 
        model.update(); 
        model.optimize(); 
        
        // Check if model was solved.  If not, LP is infeasible, not a neighbor 
        if(model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
          continue; 
        }
        
        // If model was solved, extract the constraint 
        //GRBConstr constr = model.getConstrByName("c"+std::to_string(j));
        GRBConstr constr = model.getConstrByName("c"+std::to_string(jidxs[j]));
        double constr_val = constr.get(GRB_DoubleAttr_Pi);
        
        // Return with 1 if this is a Delaunay neighbor
        if(constr_val > 0) {
          //nhblist.push_back(j); 
          nhblist.push_back(jidxs[j]); 
        }
        
        // Reset the model to unsolved state 
        model.reset(); 
      } // close j loop 
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
    
    return arma::conv_to<arma::uvec>::from(nhblist); 
  };
  

  // ***** Compute the full Delaunay adjacency, row by row *****
  struct Delaunay_ADJ_prlwkr : RcppParallel::Worker {
    // Inputs
    const GRBEnv* env;
    const arma::mat& W;
    const arma::vec& lb, ub;
    
    // Calculate parameters
    unsigned int nW;
    
    // Output
    arma::umat DADJ;
    
    // Constructor
    Delaunay_ADJ_prlwkr(const GRBEnv* env_, const arma::mat& W_, const arma::vec& lb, const arma::vec& ub) :
      env(env_), W(W_), lb(lb), ub(ub) {
      env = env_;
      nW = W.n_rows;
      DADJ = arma::zeros<arma::umat>(nW, nW);
    };
    

    // Test for Delaunay neighbors of j
    void neighbors_of_j(unsigned int j) {
      
      arma::uvec nhblist = Delaunay_nhbs_of_i(env, W, j, lb, ub, false); 
      
      for(unsigned int k=0; k<nhblist.n_elem; ++k) {
        DADJ(j,nhblist[k]) = 1; 
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
      
      Rcpp::Rcout << "Calculating Delaunay graph edges ..." << std::endl;
      for(unsigned int j=0; j<nW; ++j) {
        Rcpp::Rcout << j+1; 
        
        neighbors_of_j(j); 
        
        if(j == (nW-1) || (j+1) % 8 == 0) {
          Rcpp::Rcout << std::endl;
        } else {
          Rcpp::Rcout << "\t";
        }
        
      }
      symmetrize_DADJ(); 
      
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
      Rcpp::Rcout << "Calculating Delaunay graph edges ..." << std::endl; 
      
      while(ncompleted < nprocess) {
        Rcpp::checkUserInterrupt(); 
        
        unsigned int from = chunk_counter * chunk_size; 
        unsigned int to = from + chunk_size; 
        to = std::min(to, nprocess); 
        
        RcppParallel::parallelFor(from, to, *this); 
        
        ncompleted += to - from; 
        
        chunk_counter++; 
        
        Rcpp::Rcout << int(double(ncompleted)/double(nprocess)*100) << "%";
        if(ncompleted == nprocess || chunk_counter % 8 == 0) {
          Rcpp::Rcout << std::endl;
        } else {
          Rcpp::Rcout << "\t";
        }
      }
      
      symmetrize_DADJ(); 
      
      // log end time, print execution time 
      auto stoptime = std::chrono::high_resolution_clock::now();
      auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
      Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
    }
    
    // Symmetrize the upper-triangular matrix, after it is computed 
    void symmetrize_DADJ() {
      DADJ += DADJ.t(); 
    }
    
  };
  
  
  // ***** Compute the full Delaunay adjacency, using the Gabriel Graph as a starting seed  *****
  struct Delaunay_ADJ_Gabriel_seed_prlwkr : RcppParallel::Worker {
    
    // Inputs
    const GRBEnv* env;
    const arma::mat& W;
    const arma::vec& lb, ub;
    
    // Calculate parameters
    unsigned int nW;
    arma::umat geoDIST; // stores the geodesic distance between protos in the current DADJ 
    arma::umat checked; // stores which adjacencies have been checked after each round 
    bool is_seeded; 
    
    // Output
    arma::umat DADJ;
    
    
    // Constructor
    Delaunay_ADJ_Gabriel_seed_prlwkr(const GRBEnv* env_, const arma::mat& W_, const arma::vec& lb, const arma::vec& ub) :
      env(env_), W(W_), lb(lb), ub(ub) {
      env = env_;
      nW = W.n_rows;
      
      // Initialize DADJ as empty
      DADJ = arma::zeros<arma::umat>(nW, nW);
      
      // Initialize checked to only check upper triangle 
      // 0 means not yet checked, 1 means checked (either DADJ found or not found)
      checked = arma::zeros<arma::umat>(nW, nW); 
      checked.elem(arma::trimatl_ind( arma::size(DADJ), 0 )).ones(); 
      
      is_seeded = false; 
    };
    
    
    // Seed the DADJ with the Gabriel adjacency before computing  
    void seed_DADJ(const arma::umat& seed) {
      
      arma::uvec skip_these = arma::find(seed); 
      
      // Set the DADJ entries corresponding to the non-zero entries of the input adjacency = 1
      DADJ.elem(skip_these).ones(); 
      
      // Set the DADJ entries corresponding to the non-zero entries of the input adjacency = 1
      checked.elem(skip_these).ones(); 
      
      is_seeded = true; 
    }
    
    
    // Return a list of the (i,j) indices to be checked on next round 
    arma::umat check_these_next() {
      arma::umat out = arma::ind2sub(arma::size(DADJ), arma::find(checked==0 && geoDIST==2)).t(); 
      return out; 
    }
    
    
    // Test for Delaunay neighbors of prototype j
    void neighbors_of_j(unsigned int j) {
      
      // Find the next set of (j,k) pairs to test 
      arma::umat ck = check_these_next(); 
      ck = ck.rows(arma::find(ck.col(0)==j));
      if(ck.n_rows == 0) return; 
      
      // Return a list of Delaunay neighbors of j 
      arma::uvec nhblist = Delaunay_nhbs_of_i(env, W, j, lb, ub, false, ck.col(1)); 
      
      // Populate the DADJ of these neighbors, and mark them as checked 
      for(unsigned int k=0; k<nhblist.n_elem; ++k) {
        DADJ(j,nhblist[k]) = 1; 
        DADJ(nhblist[k],j) = 1; 
      }
      
      checked.submat(ck.col(0), ck.col(1)).ones(); 
    }
    
    
    // Parallel operator
    void operator()(std::size_t begin, std::size_t end) {
      // Loop over neighbors of prototype iidx
      for(int j = int(begin); j < int(end); j++) {
        neighbors_of_j(j); 
      } // close loop over neighbors
    } // close worker
    
    // Calc a single round 
    void calc_single_round(bool parallel) {
      
      // Bust the set of indices in i_list up into chunks so we can report
      unsigned int nprocess = nW; 
      
      unsigned int chunk_size = nprocess / 4;
      chunk_size = std::max(chunk_size, VOR_PARALLEL::get_NumThreads_RcppParallel());
      chunk_size = std::min(chunk_size, nprocess);
      unsigned int chunk_counter = 0;
      unsigned int ncompleted = 0; 
      
      
      while(ncompleted < nprocess) {
        Rcpp::checkUserInterrupt();
        
        unsigned int from = chunk_counter * chunk_size;
        unsigned int to = from + chunk_size;
        to = std::min(to, nprocess);
        
        if(parallel) {
          RcppParallel::parallelFor(from, to, *this);  
        } else {
          for(unsigned int jj=from; jj<to; ++jj) {
            neighbors_of_j(jj);
          }
        }
        
        
        ncompleted += to - from;
        
        chunk_counter++;
        
        Rcpp::Rcout << int(double(ncompleted)/double(nprocess)*100) << "%";
        
        if(chunk_counter % 8 == 0) {
          Rcpp::Rcout << std::endl;
        } else {
          Rcpp::Rcout << "\t";
        }
      }
      
    }
    
    
    // Main function to call to compute full DADJ 
    void calc_all_rounds(bool parallel) {
      
      if(!is_seeded) Rcpp::stop("Must call .seed_DADJ first");
      
      // log start time 
      auto starttime = std::chrono::high_resolution_clock::now(); 
      
      
      Rcpp::Rcout << "Testing Delaunay graph edges ..." << std::endl; 
      
      // Compute the geoDIST of the seeded DADJ and determine which proto pairs to check first 
      geoDIST = SOM_UTILS_DIST::geodesicdist(DADJ, false, false);  
      arma::umat to_check = check_these_next();
      
      unsigned int cur_round = 1; 
      unsigned int num_checked = 0; 
      
      // Calculation will be done in rounds
      // Continue with the rounds until there are no more to check 
      while(to_check.n_rows > 0) {
        Rcpp::Rcout << "Round " << cur_round << " (" << to_check.n_rows << " edges)\t"; 
        
        // Perform the calculation for this round and re-compute the geoDIST of the resulting graph 
        this->calc_single_round(parallel);
        geoDIST = SOM_UTILS_DIST::geodesicdist(DADJ, false, false);  
        num_checked += to_check.n_rows; 
        
        // Update the check list 
        to_check = check_these_next(); 
        
        Rcpp::Rcout << std::endl; 
        cur_round++; 
      }
      
      Rcpp::Rcout << "Tested " << num_checked << " of " << nW*(nW-1)/2 << " possible edges" << std::endl; 
      
      // log end time, print execution time 
      auto stoptime = std::chrono::high_resolution_clock::now();
      auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
      Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
      
    }
    
  };
  
  
  
}

#endif



/*
// ***** Old version of DADJ Worker *****
// This tested each (i,j) prototype pair separately, but in parallel. 
// The new version above tests all neighbors of i together to save time 
// (only need to build the polytope defining Voronoi cell i once)
struct Delaunay_ADJ_prlwkr : RcppParallel::Worker {
  // Inputs
  GRBEnv* env;
  const arma::mat& W;
  //int iidx;
  const arma::vec& lb, ub;
  
  // Calculate parameters
  unsigned int nW;
  arma::umat skiplist; // matrix of (i,j) indices to NOT test, 0-indexed 
  bool skipsome; 
  
  // Output
  arma::umat DADJ;
  
  // Constructor
  Delaunay_ADJ_prlwkr(GRBEnv* env_, const arma::mat& W_, const arma::vec& lb, const arma::vec& ub) :
    W(W_), lb(lb), ub(ub) {
    env = env_;
    nW = W.n_rows;
    DADJ = arma::zeros<arma::umat>(nW, nW);
    skipsome = false; 
  };
  
  // Set a list of indices to skip testing, if desired   
  void set_skiplist(const arma::umat& skiplist_) {
    skiplist = skiplist_; 
    skipsome = true; 
  }
  
  
  // Test for Delaunay neighbors of j
  void neighbors_of_j(unsigned int j) {
    
    for(unsigned int k=j+1; k<nW; ++k) {
      if(skipsome) {
        // Check if the (j,k) combo is in the skiplist 
        arma::uvec isskipped = arma::find(skiplist.col(0) == j && skiplist.col(1) == k); 
        if(isskipped.n_elem > 0) continue; 
      }
      
      DADJ(j,k) = is_Delaunay_nhb(env, W, j, k, lb, ub, false);
      //DADJ[j,k] = is_Delaunay_nhb(W, j, k, lb, ub, false); 
    }
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    // Loop over neighbors of prototype iidx
    for(int j = int(begin); j < int(end); j++) {
      neighbors_of_j(j); 
    } // close loop over neighbors
  } // close worker
  
  void calc_parallel() {
    RcppParallel::parallelFor(0, nW, *this); 
    symmetrize_DADJ(); 
  }
  
  void calc_serial() {
    // log start time 
    auto starttime = std::chrono::high_resolution_clock::now(); 
    
    Rcpp::Rcout << "Calculating Delaunay graph edges ..." << std::endl;
    for(unsigned int j=0; j<nW; ++j) {
      Rcpp::Rcout << j+1; 
      
      neighbors_of_j(j); 
      
      if(j == (nW-1) || (j+1) % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
      
    }
    symmetrize_DADJ(); 
    
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
    Rcpp::Rcout << "Calculating Delaunay graph edges ..." << std::endl; 
    
    while(ncompleted < nprocess) {
      Rcpp::checkUserInterrupt(); 
      
      unsigned int from = chunk_counter * chunk_size; 
      unsigned int to = from + chunk_size; 
      to = std::min(to, nprocess); 
      
      RcppParallel::parallelFor(from, to, *this); 
      
      ncompleted += to - from; 
      
      chunk_counter++; 
      
      Rcpp::Rcout << int(double(ncompleted)/double(nprocess)*100) << "%";
      if(ncompleted == nprocess || chunk_counter % 8 == 0) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t";
      }
    }
    
    symmetrize_DADJ(); 
    
    // log end time, print execution time 
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto exectime = std::chrono::duration<float, std::ratio<60>>(stoptime - starttime); 
    Rcpp::Rcout << "Computation time: " << exectime.count() << " minutes" << std::endl;
  }
  
  // Symmetrize the upper-triangular matrix, after it is computed 
  void symmetrize_DADJ() {
    DADJ += DADJ.t(); 
  }
  
};

struct Delaunay_presolve_worker : RcppParallel::Worker {
  // Inputs
  const arma::mat& W;
  
  // Calculate parameters
  unsigned int nW;
  
  // Output
  std::vector<arma::umat> BMUlist; 
  //arma::umat DADJ;
  
  // Constructor
  Delaunay_presolve_worker(const arma::mat& W_) :
    W(W_) {
    nW = W.n_rows; 
    //DADJ = arma::zeros<arma::umat>(nW, nW);
    BMUlist.resize(nW); 
  };
  
  // Neighbors of Wj
  void neighbors_of_j(unsigned int j) {
    
    BMUlist[j].set_size(nW-(j+1), 2); 
    
    unsigned int rowcounter = 0; 
    for(unsigned int k=j+1; k<nW; ++k) {
      
      // Find midpoint between Wj and Wk 
      arma::rowvec mid = (W.row(j) + W.row(k)) / 2.0; 
      
      // Find its BMU 
      arma::mat W_ = W; 
      W_.each_row() -= mid; 
      arma::vec dist = arma::sum(arma::square(W_), 1); 
      
      arma::uvec BMUs = arma::sort_index(dist, "ascend"); 
      //DADJ(BMUs[0],BMUs[1]) = 1; 
      BMUlist[j](rowcounter,0) = BMUs[0]; 
      BMUlist[j](rowcounter,1) = BMUs[1]; 
      rowcounter++; 
    }
    
  }
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    // Loop over neighbors of prototype iidx
    for(int j = int(begin); j < int(end); j++) {
      neighbors_of_j(j); 
    } // close loop over neighbors
  } // close worker
  
  // Parallel call 
  void calc_parallel() {
    RcppParallel::parallelFor(0, nW, *this); 
  }
  
  // Serial call 
  void calc_serial() {
    for(unsigned int j=0; j<nW; ++j) {
      neighbors_of_j(j); 
    }
  }
  
  // Build the DADJ from the BMUlist 
  arma::umat build_DADJ() {
    arma::umat DADJ(nW, nW); 
    DADJ.zeros(); 
    
    for(unsigned int j=0; j<nW; ++j) {
      arma::uvec indices = arma::sub2ind( arma::size(DADJ), BMUlist[j].t() ); 
      DADJ.elem(indices).ones(); 
    }
    
    // Symmetrize 
    DADJ += DADJ.t(); 
    return DADJ; 
    
  }
  
};
*/
