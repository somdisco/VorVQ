#ifndef VORVQ_UTILS_POLYDEF_HPP
#define VORVQ_UTILS_POLYDEF_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include "vor_utils_gurobi.hpp"


namespace VOR_POLYDEF {



/*
 Function to identify redundant constraints in a polytope definition Ax <= b
 If logical_list = false, returns a vector of indices (rows of A) which are redundant
 If logical_list = true, returns a vector (length = nrow(A)) of 1/0, 1 = constraint is redundant
 */

inline std::vector<int> cpp_polytope_redundancy(const GRBEnv* env, const arma::mat& A, const arma::vec& b, bool verbose, bool logical_list) {
  
  // Get the constraints & RHS of the linear program
  int M = A.n_rows, N = A.n_cols;
  
  // Initialize output containers
  std::vector<int> remove_set;
  
  try{
    
    // Initialize a Gurobi model
    GRBModel model = GRBModel(*env);
    if(!verbose) model.getEnv().set(GRB_IntParam_OutputFlag, 0); // Turn off
    std::vector<GRBVar> x(N);
    
    // Add the variables, set all objective coefficients to 0 to start
    for(int n=0; n<N; ++n) {
      x[n] = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0.0, GRB_CONTINUOUS, "x"+std::to_string(n)); //
    }
    
    // The objective is to maximize (see Agrell paper, we are forming his dual problem)
    model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    
    // Add the constraints
    std::vector<GRBConstr> constr(M);
    for(int m=0; m<M; ++m) {
      GRBLinExpr ntot = 0;
      for(int n=0; n<N; ++n) {
        ntot += A(m,n) * x[n];
      }
      constr[m] = model.addConstr(ntot, GRB_LESS_EQUAL, b(m), "c"+std::to_string(m));
    }
    
    // Loop over each constraint, in reverse order to preserve indexing
    for(int m=M-1; m>=0; --m) {
      
      // Loop over each variable to change its coefficient
      for(int n=0; n<N; ++n) {
        // Set the objective coefficient to A(m) for this variable
        x[n].set(GRB_DoubleAttr_Obj, A(m,n));
      }
      
      // Set the RHS to += 1
      constr[m].set(GRB_DoubleAttr_RHS, b(m)+1);
      
      // Solve the model
      model.optimize();
      
      // Check if model wasn't solved.  If not, continue
      bool constr_redundant = false;
      if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        
        // Store the optimal objective val
        double optimal_val = model.get(GRB_DoubleAttr_ObjVal);
        
        // If it is NOT > that existing RHS, the constraint is redundant
        if(!(optimal_val  > b(m))) {
          constr_redundant = true;
        }
      }
      
      if(constr_redundant) {
        // Mark it for deletion
        remove_set.push_back(m);
        
        // Remove it from model
        model.remove(constr[m]);
      } else {
        // Reset the RHS back to b
        constr[m].set(GRB_DoubleAttr_RHS, b(m));
      }
      
      
    } // close loop over constraints
    
  }
  catch (GRBException e)
  {
    if(verbose) {
      Rcpp::Rcout << "Error code = " << e.getErrorCode() << std::endl;
      Rcpp::Rcout << e.getMessage() << std::endl;
    }
    
    std::vector<int> tmp(0);
    return tmp;
  }
  catch (...)
  {
    if(verbose) {
      Rcpp::Rcout << "Exception during optimization" << std::endl;
    }
    
    std::vector<int> tmp(0);
    return tmp;
  }
  
  // Determine what type of return we want 
  if(!logical_list) {
    return remove_set;  
  } else {
    
    // Initialize a vector of length = # of constraints, every entry initialized = 0 
    std::vector<int> constr_is_redundant(M);
    std::fill(constr_is_redundant.begin(), constr_is_redundant.end(), 0);
    // Put a 1 in the indices to be removed
    for(unsigned int k=0; k<remove_set.size(); ++k) {
      constr_is_redundant[remove_set[k]] = 1; 
    }
    
    return constr_is_redundant; 
  } 
}

/*
 Function to drop redundant rows for a polytope system 
 Several overloaded versions
 */
inline void cpp_polytope_rmv_redundancy(const arma::mat& A_in, const arma::vec& b_in, arma::mat& A_out, arma::vec& b_out) {
  
  // Initialize a Gurobi environment
  // FILE* mystdout = stdout;
  // stdout = fopen("/dev/null", "w");
  // std::fclose(stdout);
  // GRBEnv* env = 0;
  // env = new GRBEnv();
  // stdout = mystdout;
  // 
  // Get list of logical redundancies (vector of length = #constraints, 1=redundant, 0=not)
  // and calculate how many constraints are redundant 
  std::vector<int> constr_is_redundant = cpp_polytope_redundancy(GUROBIENV, A_in, b_in, false, true); // verbose=F, logical=T
  int nredundant = std::accumulate(constr_is_redundant.begin(), constr_is_redundant.end(), 0);
  
  if(nredundant > 0) {
    // Keep the constraints that have a logical flag = 0 
    arma::uvec keep_these(A_in.n_rows - nredundant); 
    int counter = 0; 
    for(unsigned int i=0; i<A_in.n_rows; ++i) {
      if(constr_is_redundant[i] > 0) continue;
      keep_these[counter] = i;
      counter++;
    }
    
    // Copy the keep rows 
    A_out = A_in.rows(keep_these);
    b_out = b_in.rows(keep_these);
  } else {
    // Just set the out container = in containers 
    A_out = A_in;
    b_out = b_in;
  }
  
  // Clean up Gurobi environment
  //delete env;
}

inline void cpp_polytope_rmv_redundancy(const arma::mat& A_in, const arma::vec& b_in, const arma::uvec& cid_in, arma::mat& A_out, arma::vec& b_out, arma::uvec& cid_out) {
  
  // Initialize a Gurobi environment
  // FILE* mystdout = stdout;
  // stdout = fopen("/dev/null", "w");
  // std::fclose(stdout);
  // GRBEnv* env = 0;
  // env = new GRBEnv();
  // stdout = mystdout;
  
  // Get list of logical redundancies (vector of length = #constraints, 1=redundant, 0=not)
  // and calculate how many constraints are redundant 
  std::vector<int> constr_is_redundant = cpp_polytope_redundancy(GUROBIENV, A_in, b_in, false, true); // verbose=F, logical=T
  int nredundant = std::accumulate(constr_is_redundant.begin(), constr_is_redundant.end(), 0);
  
  if(nredundant > 0) {
    // Keep the constraints that have a logical flag = 0 
    arma::uvec keep_these(A_in.n_rows - nredundant); 
    int counter = 0; 
    for(unsigned int i=0; i<A_in.n_rows; ++i) {
      if(constr_is_redundant[i] > 0) continue;
      keep_these[counter] = i;
      counter++;
    }
    
    // Copy the keep rows 
    A_out = A_in.rows(keep_these);
    b_out = b_in.rows(keep_these);
    cid_out = cid_in.rows(keep_these);
  } else {
    // Just set the out container = in containers 
    A_out = A_in;
    b_out = b_in;
    cid_out = cid_in;
  }
  
  // Clean up Gurobi environment
  //delete env;
}






// H-representation of First-Order Voronoi Cell,
// WITHOUT global bounds appended to the constraints
// A, b, cid are output containers
inline void cpp_vor1_polytope(const arma::mat& W, int iidx, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  int nW = W.n_rows, d = W.n_cols;
  A.resize(nW-1, d);
  b.resize(nW-1);
  cid.resize(nW-1);
  
  double norm2_iidx = arma::accu(arma::square(W.row(iidx)));
  
  int counter = 0;
  for(int k=0; k<nW; ++k) {
    if(k == iidx) continue;
    A.row(counter) = W.row(k) - W.row(iidx);
    
    double norm2_k = arma::accu(arma::square(W.row(k)));
    b(counter) = (norm2_k - norm2_iidx) / 2.0;
    
    cid(counter) = k;
    counter++;
  }
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}

// H-representation of First-Order Voronoi cell,
// WITH global x bounds appended to the contraints
// A, b, cid are output containers
inline void cpp_vor1_polytope(const arma::mat& W, int iidx, const arma::vec& xlb, const arma::vec& xub, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  
  // First set A, b, cid
  cpp_vor1_polytope(W, iidx, A, b, cid);
  
  // Now we need to add two diagonal (d x d) matrices to the end of A:
  // One for the global upper bound (Ax <= xub) and one for the global lower bound (-Ax <= -xlb)
  int nW = W.n_rows, d = W.n_cols;
  arma::mat A_append(d,d); A_append.eye();
  arma::uvec cid_append(d); cid_append.fill(nW);
  
  // Add the global upper bound
  A = arma::join_cols(A, A_append);
  b = arma::join_cols(b, xub);
  cid = arma::join_cols(cid, cid_append);
  
  // Add the global lower bound
  A = arma::join_cols(A, -A_append);
  b = arma::join_cols(b, -xlb);
  cid = arma::join_cols(cid, cid_append);
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}

// H-Representation of First-Order Voronoi Cell,
// where redundant constraints are removed using the Delaunay Adjacency,
// WITHOUT global bounds appended to the constraints
inline void cpp_vor1_polytope(const arma::mat& W, const arma::umat& DADJ, int iidx, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  
  // Determine the Delaunay neighbors of i
  arma::uvec nhbs_of_i = arma::find(DADJ.row(iidx));
  if(arma::any(nhbs_of_i==iidx)) { // Make sure that i is not its own neighbor
    nhbs_of_i = nhbs_of_i.elem(arma::find(nhbs_of_i != iidx));
  }
  int nnhb = nhbs_of_i.size();
  
  // Initialize output containers to hold the 1st-order & 2nd-order cell constraints
  A.resize(nnhb, W.n_cols);
  b.resize(nnhb);
  cid.resize(nnhb);
  
  // Precompute squared norms of i & j prototypes
  double norm2_iidx = arma::accu(arma::square(W.row(iidx)));
  
  // Loop
  int counter = 0;
  for(int k=0; k<nnhb; ++k) {
    if(nhbs_of_i(k) == unsigned(iidx)) continue;
    
    A.row(counter) = W.row(nhbs_of_i(k)) - W.row(iidx);
    double norm2_k = arma::accu(arma::square( W.row(nhbs_of_i(k)) ));
    b(counter) = (norm2_k - norm2_iidx) / 2.0;
    
    cid(counter) = nhbs_of_i(k);
    
    counter++;
  }
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}

// H-Representation of First-Order Voronoi Cell,
// where redundant constraints are removed using the Delaunay Adjacency,
// WITH global bounds appended to the constraints
inline void cpp_vor1_polytope(const arma::mat& W, const arma::umat& DADJ, int iidx, const arma::vec& xlb, const arma::vec& xub, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  
  // First set A, b, cid
  cpp_vor1_polytope(W, DADJ, iidx, A, b, cid);
  
  // Now we need to add two diagonal (d x d) matrices to the end of A:
  // One for the global upper bound (Ax <= xub) and one for the global lower bound (-Ax <= -xlb)
  int nW = W.n_rows, d = W.n_cols;
  arma::mat A_append(d,d); A_append.eye();
  arma::uvec cid_append(d); cid_append.fill(nW);
  
  // Add the global upper bound
  A = arma::join_cols(A, A_append);
  b = arma::join_cols(b, xub);
  cid = arma::join_cols(cid, cid_append);
  
  // Add the global lower bound
  A = arma::join_cols(A, -A_append);
  b = arma::join_cols(b, -xlb);
  cid = arma::join_cols(cid, cid_append);
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}




// H-Representation of Second-Order Voronoi Cell,
// where redundant constraints are removed using the Delaunay Adjacency
// WITHOUT global bounds appended to the constraints
inline void cpp_vor2_polytope(const arma::mat& W, const arma::umat& DADJ, int iidx, int jidx, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  
  // ** Check that iidx & jidx are Delaunay neighbors. If not, can't proceed.
  if(!(DADJ(iidx,jidx) > 0)) Rcpp::stop("Prototypes iidx & jidx are not Delaunay adjacent.");
  
  // Determine the Delaunay neighbors of i
  arma::uvec nhbs_of_i = arma::find(DADJ.row(iidx));
  if(arma::any(nhbs_of_i==iidx)) { // Make sure that i is not its own neighbor
    nhbs_of_i = nhbs_of_i.elem(arma::find(nhbs_of_i != iidx));
  }
  int nnhb = nhbs_of_i.size();
  
  // Initialize output containers to hold the 1st-order & 2nd-order cell constraints
  A.resize(nnhb + nnhb-1, W.n_cols);
  b.resize(nnhb + nnhb-1);
  cid.resize(nnhb + nnhb-1);
  
  // Precompute squared norms of i & j prototypes
  double norm2_iidx = arma::accu(arma::square(W.row(iidx)));
  double norm2_jidx = arma::accu(arma::square(W.row(jidx)));
  
  // Loop
  int counter = 0;
  for(int k=0; k<nnhb; ++k) {
    if(nhbs_of_i(k) == unsigned(iidx)) continue;
    
    A.row(counter) = W.row(nhbs_of_i(k)) - W.row(iidx);
    double norm2_k = arma::accu(arma::square( W.row(nhbs_of_i(k)) ));
    b(counter) = (norm2_k - norm2_iidx) / 2.0;
    
    cid(counter) = nhbs_of_i(k);
    
    counter++;
    
    if(nhbs_of_i(k) == unsigned(jidx)) continue;
    
    A.row(counter) = W.row(nhbs_of_i(k)) - W.row(jidx);
    b(counter) = (norm2_k - norm2_jidx) / 2.0;
    
    cid(counter) = nhbs_of_i(k);
    
    counter++;
  }
  
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}

// H-Representation of Second-Order Voronoi Cell,
// where redundant constraints are removed using the Delaunay Adjacency
// WITHOUT global bounds appended to the constraints
inline void cpp_vor2_polytope(const arma::mat& W, const arma::umat& DADJ, int iidx, int jidx, const arma::vec& xlb, const arma::vec& xub, arma::mat& A, arma::vec& b, arma::uvec& cid, bool rmv_redundancies = false) {
  
  // ** Check that iidx & jidx are Delaunay neighbors. If not, can't proceed.
  if(!(DADJ(iidx,jidx) > 0)) Rcpp::stop("Prototypes iidx & jidx are not Delaunay adjacent.");
  
  // First set A, b, cid
  cpp_vor2_polytope(W, DADJ, iidx, jidx, A, b, cid);
  
  // Now we need to add two diagonal (d x d) matrices to the end of A:
  // One for the global upper bound (Ax <= xub) and one for the global lower bound (-Ax <= -xlb)
  int nW = W.n_rows, d = W.n_cols;
  arma::mat A_append(d,d); A_append.eye();
  arma::uvec cid_append(d); cid_append.fill(nW);
  
  // Add the global upper bound
  A = arma::join_cols(A, A_append);
  b = arma::join_cols(b, xub);
  cid = arma::join_cols(cid, cid_append);
  
  // Add the global lower bound
  A = arma::join_cols(A, -A_append);
  b = arma::join_cols(b, -xlb);
  cid = arma::join_cols(cid, cid_append);
  
  
  if(rmv_redundancies) {
    cpp_polytope_rmv_redundancy(A, b, cid, A, b, cid);
  }
}




} // close namespace

#endif