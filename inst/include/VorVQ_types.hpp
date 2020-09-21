#ifndef VORVQ_VOROBJ_HPP
#define VORVQ_VOROBJ_HPP


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <algorithm>
//#include <typeinfo>


#include "vor_utils_gabriel.hpp"
#include "vor_utils_delaunay.hpp"
#include "vor_utils_mvie.hpp"
#include "vor_utils_polyinterior.hpp"
#include "vor_utils_dikin.hpp"


// Forward declaration for RCPP_MODULES 
class VOR; 
RCPP_EXPOSED_CLASS(VOR);


// Container for MVIE parameters 
struct VOR_params_container {
  
  bool parallel; 
  
  int MVIE_maxiter; 
  double MVIE_tol; 
  bool MVIE_fix_c0; 
  
  bool seed_DADJ; 
  
  // Set default values in constructor 
  VOR_params_container() {
    parallel = true; 
    
    MVIE_maxiter = 100; 
    MVIE_tol = 1e-4; 
    MVIE_fix_c0 = false; 
    
    seed_DADJ = true; 
  }
  
  void set(Rcpp::List paramlist) {
    
    Rcpp::CharacterVector paramname = paramlist.names(); 
    
    for(unsigned int i=0; i<paramlist.size(); ++i) {
      if(paramname[i] == "parallel") {
        this->parallel = Rcpp::as<bool>(paramlist[i]); 
        continue;
      }
      
      // if(paramname[i] == "polydef_ADJ") {
      //   std::string tmp = Rcpp::as<std::string>(paramlist[i]); 
      //   if(tmp=="none" || tmp=="dadj" || tmp=="cadj") {
      //     this->polydef_ADJ = tmp; 
      //   } else {
      //     Rcpp::stop("polydef_ADJ must be one of 'none','dadj','cadj'");
      //   }
      //   continue; 
      // }
      
      if(paramname[i] == "MVIE_maxiter") {
        int tmp = Rcpp::as<int>(paramlist[i]); 
        if(tmp > 0) {
          this->MVIE_maxiter = tmp; 
        } else {
          Rcpp::stop("MVIE_maxiter must be > 0");
        }
        continue; 
      }
      
      if(paramname[i] == "MVIE_tol") {
        double tmp = Rcpp::as<int>(paramlist[i]); 
        if(tmp > 0) {
          this->MVIE_tol = tmp; 
        } else {
          Rcpp::stop("MVIE_tol must be > 0");
        }
        continue; 
      }
      
      if(paramname[i] == "MVIE_fix_c0") {
        this->MVIE_fix_c0 = Rcpp::as<bool>(paramlist[i]); 
        continue;
      }
      
      // if(paramname[i] == "vor1_centers") {
      //   std::string tmp = Rcpp::as<std::string>(paramlist[i]); 
      //   if(tmp=="W" || tmp=="interior" || tmp=="chebyshev" || tmp=="analytic") {
      //     this->vor1_centers = tmp; 
      //   } else {
      //     Rcpp::warning("Unknown value for vor1_centers. Options are 'W', 'interior', 'chebyshev', 'analytic'");
      //   }
      //   continue; 
      // }
      
      // if(paramname[i] == "vor2_centers") {
      //   std::string tmp = Rcpp::as<std::string>(paramlist[i]); 
      //   if(tmp=="interior" || tmp=="chebyshev" || tmp=="analytic") {
      //     this->vor1_centers = tmp; 
      //   } else {
      //     Rcpp::warning("Unknown value for vor1_centers. Options are 'interior', 'chebyshev', 'analytic'");
      //   }
      //   continue; 
      // }
      
      if(paramname[i] == "seed_DADJ") {
        this->seed_DADJ = Rcpp::as<bool>(paramlist[i]); 
        continue;
      }
      
      std::string warnstr = "Unknown parameter: " + paramname[i]; 
      warnstr += "\nShould be one of 'parallel','polydef_ADJ','MVIE_maxiter','MVIE_tol','MVIE_fix_c0','vor1_centers','vor2_centers'";
      Rcpp::warning(warnstr); 
    }

    
    return; 
    
  } // close set 
  
  
  Rcpp::List get() {
    Rcpp::List out;
    out["parallel"] = parallel; 
    //out["polydef_ADJ"] = polydef_ADJ; 
    out["MVIE_maxiter"] = MVIE_maxiter; 
    out["MVIE_tol"] = MVIE_tol; 
    out["MVIE_fix_c0"] = MVIE_fix_c0; 
    //out["vor1_centers"] = vor1_centers; 
    //out["vor2_centers"] = vor2_centers; 
    out["seed_DADJ"] = seed_DADJ; 
    return out; 
  }
  
};



class VOR{
  
public:
  
  
  // Constructor
  VOR(); 
  
  // Control parameters 
  VOR_params_container params; 
  void set_params(Rcpp::List paramlist); 
  Rcpp::List get_params(); 

  
  // Prototypes & bounds 
  arma::mat W; 
  unsigned int nW; 
  unsigned int d; 
  arma::vec lb; 
  arma::vec ub; 
  arma::umat CADJ; 
  arma::umat GADJ; 
  
  void set_W(const arma::mat& W_); 
  void set_bounds(arma::vec lb_, arma::vec ub_); 
  void set_CADJ(const arma::umat& CADJ_); 
  void calc_GADJ(); 
  void initialize_VOR(const arma::mat& W_, const arma::umat& CADJ_); 
  
  // List of active vor1 and vor2 indices to use 
  // Stored as 1-based 
  arma::uvec vor1_active; 
  arma::umat vor2_active; 
  
  void set_vor1_active(const arma::uvec& vor1); 
  void set_vor2_active(const arma::umat& vor2);
  
  // Polytope defs 
  Rcpp::List get_vor1_polytope(unsigned int iidx); 
  Rcpp::List get_vor2_polytope(unsigned int iidx, unsigned int jidx); 
  
  // Delauanay
  arma::umat DADJ; 
  
  void calc_DADJ(); 
  void set_DADJ(const arma::umat& DADJ_); 
  
  
  // Chebyshev centers 
  arma::mat vor1_centers; 
  arma::mat vor2_centers;

  void calc_vor1_centers(); 
  void calc_vor2_centers(); 
  void set_vor1_centers(const arma::mat& C); 
  void set_vor2_centers(const arma::mat& C); 
  arma::mat get_vor1_centers(); 
  arma::mat get_vor2_centers(); 
  void clear_vor1_centers(); 
  void clear_vor2_centers(); 

  // MVIE 
  arma::mat vor1_MVIE_c; 
  arma::cube vor1_MVIE_E; 
  arma::vec vor1_MVIE_logdet;
  arma::vec vor1_MVIE_volratio; 
  arma::ivec vor1_MVIE_status;
  
  arma::mat vor2_MVIE_c; 
  arma::cube vor2_MVIE_E; 
  arma::vec vor2_MVIE_logdet; 
  arma::vec vor2_MVIE_volratio; 
  arma::ivec vor2_MVIE_status; 
  
  void calc_vor1_MVIE(); 
  void calc_vor2_MVIE(); 
  
  
  // Dikin 
  arma::cube vor1_Dikin_E; 
  arma::cube vor2_Dikin_E; 
  
  void calc_vor1_Dikin(); 
  void calc_vor2_Dikin(); 
  
  void calc_all(); 
  
  // Vis 
  Rcpp::List vis_par; 
  void set_vis_par(Rcpp::List par_); 
  
  // Flags
  bool isinit; 
  bool isset_GADJ; 
  bool isset_DADJ; 
  
  bool isset_vor2_centers; 
  
  bool isset_vor1_MVIE;
  bool isset_vor2_MVIE; 
  
  bool isset_vor1_Dikin; 
  bool isset_vor2_Dikin; 
  
  // IO 
  void save(std::string rdsfile);
  void load(std::string rdsfile); 
  
};


// ***** Constructor 
inline VOR::VOR() {
  
  // Control
  this->params = VOR_params_container(); 

  
  this->vor1_active.set_size(0); 
  this->vor2_active.set_size(0,0); 
  
  this->isinit = false; 
  this->isset_GADJ = false; 
  this->isset_DADJ = false; 
  this->isset_vor1_MVIE = false; 
  this->isset_vor2_MVIE = false; 
  this->isset_vor2_centers = false; 
  this->isset_vor1_Dikin = false; 
  this->isset_vor2_Dikin = false; 
  
}

inline void VOR::set_bounds(arma::vec lb_, arma::vec ub_) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  
  // Check sizes 
  if(lb_.n_elem == 1) {
    arma::vec tmplb = arma::vec(this->d); 
    tmplb.fill(lb_[0]); 
    lb_ = tmplb; 
  } else if(lb_.n_elem != d) {
    Rcpp::stop("length(lb) must = 1 or ncol(W)");
  }
  
  if(ub_.n_elem == 1) {
    arma::vec tmpub = arma::vec(this->d); 
    tmpub.fill(ub_[0]); 
    ub_ = tmpub; 
  } else if(ub_.n_elem != d) {
    Rcpp::stop("length(ub) must = 1 or ncol(W)");
  }
  
  // Check values 
  if(arma::any(lb_ > ub_)) Rcpp::stop("All lb must be <= ub");
  
  arma::rowvec Wmin = arma::min(W, 0); 
  arma::rowvec Wmax = arma::max(W, 0); 
  if(arma::any(lb_ > Wmin.t())) Rcpp::stop("All lb must be <= min(W)");
  if(arma::any(ub_ < Wmax.t())) Rcpp::stop("All ub must be >= max(W)");
  
  this->lb = lb_; 
  this->ub = ub_; 
}

inline void VOR::calc_GADJ() {
  
  // Initialize parallel structure & call it 
  VOR_GABRIEL::Gabriel_ADJ_prlwkr wkr(W, 0); 
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  this->GADJ = wkr.GADJ; 
  this->isset_GADJ = true; 
} 

inline void VOR::initialize_VOR(const arma::mat& W_, const arma::umat& CADJ_) {
  
  // Set W 
  Rcpp::Rcout << "Initializing Voronoi object" << std::endl;
  
  Rcpp::Rcout << "++ storing W and dimensions ... ";
  this->W = W_; 
  this->nW = W_.n_rows; 
  this->d = W_.n_cols; 
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "++ setting global bounds = +/- 5% of W range ... ";
  this->lb = arma::min(this->W, 0).t(); 
  this->ub = arma::max(this->W, 0).t(); 
  arma::vec rng = this->ub - this->lb; 
  this->lb -= 0.05*rng; 
  this->ub += 0.05*rng; 
  Rcpp::Rcout << "done" << std::endl; 
  Rcpp::Rcout << "   call $set_bounds to change" << std::endl; 
  
  
  Rcpp::Rcout << "++ storing CADJ ... "; 
  if(CADJ_.n_rows != this->nW || CADJ_.n_cols != this->nW) Rcpp::stop("nrow(CADJ) or ncol(CADJ) != nW");
  this->CADJ = CADJ_; 
  Rcpp::Rcout << "done" << std::endl; 
  
  
  Rcpp::Rcout << "++ setting active vor1 indices (default to non-empty vor1 cells) ... ";
  arma::uvec CADJsumempty = arma::find(arma::sum(this->CADJ, 1)); // rowsums 
  CADJsumempty += 1; 
  this->vor1_active = CADJsumempty; 
  Rcpp::Rcout << "done" << std::endl; 
  Rcpp::Rcout << "   call $set_vor1_active to change" << std::endl; 
  
  
  Rcpp::Rcout << "++ setting active vor2 indices (default to non-empty vor2 cells) ... ";
  arma::umat CADJnonempty = arma::ind2sub(arma::size(this->CADJ), arma::find(this->CADJ)).t(); 
  CADJnonempty += 1;
  this->vor2_active = CADJnonempty; 
  Rcpp::Rcout << "done" << std::endl; 
  Rcpp::Rcout << "   call $set_vor2_active to change" << std::endl; 
  
  this->calc_GADJ(); 
  
  this->isinit = true; 
}


// ***** Control 
inline void VOR::set_params(Rcpp::List paramlist) {
  this->params.set(paramlist); 
}

inline Rcpp::List VOR::get_params() {
  return this->params.get(); 
}



// ***** Active Vor indices 
inline void VOR::set_vor1_active(const arma::uvec& vor1) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(arma::any(vor1 > this->nW)) Rcpp::stop("vor1 index list contains out-of-bound entries");
  
  this->vor1_active = vor1;
}

inline void VOR::set_vor2_active(const arma::umat& vor2) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(vor2.n_cols != 2) Rcpp::stop("ncol(vor2 index list) != 2");
  if(arma::any(arma::vectorise(vor2) > this->nW)) Rcpp::stop("vor2 index list contains out-of-bound entries");
  
  // Check that vor2 cells are Delaunay adjacency, if DADJ has been calculated 
  if(this->isset_DADJ) {
    bool stopthis = false; 
    for(unsigned int i=0; i<vor2.n_rows; ++i) {
      if(DADJ(vor2(i,0)-1, vor2(i,1)-1) == 0) {
        Rcpp::Rcout << "vor2 cell (" << vor2(i,0) << "," << vor2(i,1) << ") is not Delaunay adjacenct" << std::endl; 
        stopthis = true; 
      }
    }
    
   if(stopthis) Rcpp::stop("All input vor2 cells must be Delaunay adjacent"); 
  }
  
  this->vor2_active = vor2; 
}



// ***** Polytope defs 
inline Rcpp::List VOR::get_vor1_polytope(unsigned int iidx) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  arma::mat A; arma::vec b; arma::uvec cid; 
  VOR_POLYDEF::cpp_vor1_polytope(this->W, this->DADJ, iidx-1, this->lb, this->ub, A, b, cid, false); 
  Rcpp::List out; 
  out["A"] = A; 
  out["b"] = b; 
  return out; 
}

inline Rcpp::List VOR::get_vor2_polytope(unsigned int iidx, unsigned int jidx) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  if(!(this->DADJ(iidx-1,jidx-1)==1)) Rcpp::stop("DADJ(iidx,jidx)!=1: Voronoi cells iidx and jidx are not adjacenct, no second-order Voronoi polytope exists.");
  
  arma::mat A; arma::vec b; arma::uvec cid; 
  VOR_POLYDEF::cpp_vor2_polytope(this->W, this->DADJ, iidx-1, jidx-1, this->lb, this->ub, A, b, cid, false); 
  Rcpp::List out; 
  out["A"] = A; 
  out["b"] = b; 
  return out; 
}



// ***** Delaunay 
inline void VOR::calc_DADJ() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  
  if(params.seed_DADJ) {
    // Initialize parallel structure & call it 
    VOR_DELAUNAY::Delaunay_ADJ_Gabriel_seed_prlwkr wkr(GUROBIENV, W, lb, ub);
    wkr.seed_DADJ(this->GADJ); 
    wkr.seed_DADJ(this->CADJ); 
    wkr.seed_DADJ(this->CADJ.t());
    wkr.calc_all_rounds(params.parallel); 
    this->DADJ = wkr.DADJ; 
  } else {
    // Initialize parallel structure & call it 
    VOR_DELAUNAY::Delaunay_ADJ_prlwkr wkr(GUROBIENV, W, lb, ub);
    if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
    this->DADJ = wkr.DADJ;  
  }
  
  this->isset_DADJ = true; 
}

inline void VOR::set_DADJ(const arma::umat& DADJ_) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(DADJ_.n_rows != this->nW || DADJ_.n_cols != this->nW) Rcpp::stop("nrow(DADJ) or ncol(DADJ) != nW");
  
  this->DADJ = DADJ_; 
  this->isset_DADJ = true; 
}


// ***** Chebyshev centers 
inline void VOR::calc_vor1_centers() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
 
  // Initialize the struct
  arma::uvec i_list = this->vor1_active; 
  i_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor1_chebyshev_center_prlwkr wkr(this->W, this->DADJ, i_list, this->lb, this->ub); 
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Otherwise store the points 
  this->vor1_centers = wkr.C; 
}

inline void VOR::calc_vor2_centers() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor2_active.n_elem == 0) Rcpp::stop("Must call $set_vor2_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  // Initialize the struct
  arma::umat ij_list = this->vor2_active; 
  ij_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor2_chebyshev_center_prlwkr wkr(this->W, this->DADJ, ij_list, this->lb, this->ub);
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Store the points 
  this->vor2_centers = wkr.C; 
  this->isset_vor2_centers = true;
}


// ***** User-specified centers 
inline void VOR::set_vor1_centers(const arma::mat& C) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  // Check that C has nrows = length(vor1_active)
  if(C.n_rows != this->vor1_active.n_elem) Rcpp::stop("nrow(C) must equal length(vor1_active)");
  
  // Initialize the struct
  arma::uvec i_list = this->vor1_active; 
  i_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor1_check_interior_prlwkr wkr(this->W, this->DADJ, i_list, this->lb, this->ub, C); 
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Check if any points are not interior, if so print their indices to console and stop. 
  if(arma::any(wkr.is_interior == 0)) {
    arma::uvec not_interior = arma::find(wkr.is_interior==0); 
    not_interior += 1; 
    
    Rcpp::Rcout << "These rows of C are not interior to their vor1 polytopes:" << std::endl; 
    for(unsigned int i = 0; i < not_interior.n_elem; ++i) {
      Rcpp::Rcout << not_interior[i]; 
      if((i+1) % 8 == 0 || (i+1)==not_interior.n_elem) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t"; 
      }
    }
    
    Rcpp::stop("All points in C must be interior");
  }
  
  // Otherwise store the points 
  this->vor1_centers = C; 
}

inline arma::mat VOR::get_vor1_centers() {
  
  // Returns the values stored in vor1_centers if they have been set, 
  // otherwise returns the active rows of W 
  if(this->vor1_centers.n_rows == this->vor1_active.n_elem) {
    return this->vor1_centers; 
  } else { 
    arma::uvec i_list = this->vor1_active; 
    i_list -= 1; 
    return this->W.rows(i_list); 
  }
}

inline arma::mat VOR::get_vor2_centers() {
  
  // Returns the values stored in vor2_centers if they have been set, 
  // otherwise returns an error 
  if(this->vor2_centers.n_rows == this->vor2_active.n_rows) {
    return this->vor2_centers; 
  } else { 
    Rcpp::stop("Must call $calc/set_vor2_centers before calling $get_vor2_centers");
  }
}

inline void VOR::set_vor2_centers(const arma::mat& C) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor2_active.n_elem == 0) Rcpp::stop("Must call $set_vor2_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  // Check that C has nrows = nrow(vor2_active)
  if(C.n_rows != this->vor2_active.n_rows) Rcpp::stop("nrow(C) must equal nrow(vor2_active)");
  
  
  // Initialize the struct
  arma::umat ij_list = this->vor2_active; 
  ij_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor2_check_interior_prlwkr wkr(this->W, this->DADJ, ij_list, this->lb, this->ub, C); 
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Check if any points are not interior, if so print their indices to console and stop. 
  if(arma::any(wkr.is_interior == 0)) {
    arma::uvec not_interior = arma::find(wkr.is_interior==0); 
    not_interior += 1; 
    
    Rcpp::Rcout << "These rows of C are not interior to their vor2 polytopes:" << std::endl; 
    for(unsigned int i = 0; i < not_interior.n_elem; ++i) {
      Rcpp::Rcout << not_interior[i]; 
      if((i+1) % 8 == 0 || (i+1)==not_interior.n_elem) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t"; 
      }
    }
    
    Rcpp::stop("All points in C must be interior");
  }
  
  // Otherwise store the points 
  this->vor2_centers = C; 
}

inline void VOR::clear_vor1_centers() {
  this->vor1_centers.set_size(0,0); 
}

inline void VOR::clear_vor2_centers() {
  this->vor2_centers.set_size(0,0); 
}


// // ***** MVIE 
inline void VOR::calc_vor1_MVIE() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");

  // Initialize the struct
  arma::uvec i_list = this->vor1_active;
  i_list -= 1; // -1 since we're passing to c++
  
  arma::mat C0 = this->get_vor1_centers();
  VOR_MVIE::vor1_mvie_solve_prlwkr wkr(this->W, this->DADJ, i_list,
                                       this->lb, this->ub, C0,
                                       this->params.MVIE_maxiter, this->params.MVIE_tol, this->params.MVIE_fix_c0,
                                       true, true); // last is use_bounds flag

  // Run, according to parallel setting
  if(this->params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial();

  // Store
  this->vor1_MVIE_c = wkr.C;
  this->vor1_MVIE_E = wkr.E;
  this->vor1_MVIE_status = wkr.status;
  this->vor1_MVIE_logdet = wkr.logdet;
  this->vor1_MVIE_volratio = wkr.volratio;

  this->isset_vor1_MVIE = true;
}

inline void VOR::calc_vor2_MVIE() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor2_active.n_elem == 0) Rcpp::stop("Must call $set_vor2_active first.");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  if(!this->isset_vor2_centers) Rcpp::stop("Must call $calc_vor2_centers first.");
  
  
  // Initialize the struct
  arma::umat ij_list = this->vor2_active;
  ij_list -= 1; // -1 since we're passing to c++
  
  VOR_MVIE::vor2_mvie_solve_prlwkr wkr(this->W, this->DADJ, ij_list,
                                       this->lb, this->ub, this->vor2_centers,
                                       this->params.MVIE_maxiter, this->params.MVIE_tol, this->params.MVIE_fix_c0);

  // Run, according to parallel setting
  if(this->params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial();

  // Store
  this->vor2_MVIE_c = wkr.C;
  this->vor2_MVIE_E = wkr.E;
  this->vor2_MVIE_status = wkr.status;
  this->vor2_MVIE_logdet = wkr.logdet;
  this->vor2_MVIE_volratio = wkr.volratio;

  this->isset_vor2_MVIE = true;
}


// ***** Dikin
inline void VOR::calc_vor1_Dikin() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  // Initialize the struct
  arma::uvec i_list = this->vor1_active;
  i_list -= 1; // -1 since we're passing to c++
  
  VOR_DIKIN::vor1_dikin_ell_prlwkr wkr(this->W, this->DADJ, this->lb, this->ub, i_list, this->get_vor1_centers());

  // Run, according to parallel setting
  if(this->params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial();

  // Store
  this->vor1_Dikin_E = wkr.E;

  this->isset_vor1_Dikin = true;
}

inline void VOR::calc_vor2_Dikin() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor2_active.n_elem == 0) Rcpp::stop("Must call $set_vor2_active first.");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  if(!this->isset_vor2_centers) Rcpp::stop("Must call $calc_vor2_centers first.");
  
  // Initialize the struct
  arma::umat ij_list = this->vor2_active;
  ij_list -= 1; // -1 since we're passing to c++
  
  VOR_DIKIN::vor2_dikin_ell_prlwkr wkr(this->W, this->DADJ, this->lb, this->ub, ij_list, this->vor2_centers);

  // Run, according to parallel setting
  if(this->params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial();

  // Store
  this->vor2_Dikin_E = wkr.E;

  this->isset_vor2_Dikin = true;
}


// ***** Calc everything 
inline void VOR::calc_all() {
  // DADJ 
  this->calc_DADJ(); 
  Rcpp::Rcout << std::endl; 
  
  // Centers 
  this->calc_vor2_centers();
  Rcpp::Rcout << std::endl; 
  
  // Dikin 
  this->calc_vor1_Dikin(); 
  Rcpp::Rcout << std::endl; 
  this->calc_vor2_Dikin(); 
  Rcpp::Rcout << std::endl; 
  
  // MVIE 
  this->calc_vor1_MVIE(); 
  Rcpp::Rcout << std::endl; 
  this->calc_vor2_MVIE(); 
  Rcpp::Rcout << std::endl; 
}


// ***** Vis 
inline void VOR::set_vis_par(Rcpp::List par_) {
  this->vis_par = par_; 
}


// ***** Input / Output 
inline void VOR::save(std::string rdsfile) {
  
  std::string rdsend = ".vor"; 
  if(rdsend.size() > rdsfile.size() || !std::equal(rdsend.rbegin(), rdsend.rend(), rdsfile.rbegin()))
    Rcpp::stop("Save file string must have extension .vor");
  
  Rcpp::List VORList; 
  
  VORList["params_parallel"] = this->params.parallel; 
  VORList["params_MVIE_maxiter"] = this->params.MVIE_maxiter; 
  VORList["params_MVIE_tol"] = this->params.MVIE_tol; 
  VORList["params_MVIE_fix_c0"] = this->params.MVIE_fix_c0; 
  VORList["params_seed_DADJ"] = this->params.seed_DADJ; 
  
  // Prototypes & bounds 
  VORList["W"] = this->W; 
  VORList["nW"] = this->nW; 
  VORList["d"] = this->d; 
  VORList["lb"] = this->lb; 
  VORList["ub"] = this->ub; 
  VORList["CADJ"] = this->CADJ; 
  VORList["GADJ"] = this->GADJ; 
  
  VORList["vor1_active"] = this->vor1_active; 
  VORList["vor2_active"] = this->vor2_active; 
  
  VORList["DADJ"] = this->DADJ; 
  
  VORList["vor1_centers"] = this->vor1_centers; 
  VORList["vor2_centers"] = this->vor2_centers;
  
  VORList["vor1_MVIE_c"] = this->vor1_MVIE_c; 
  VORList["vor1_MVIE_E"] = this->vor1_MVIE_E; 
  VORList["vor1_MVIE_logdet"] = this->vor1_MVIE_logdet;
  VORList["vor1_MVIE_volratio"] = this->vor1_MVIE_volratio; 
  VORList["vor1_MVIE_status"]  = this->vor1_MVIE_status;
  
  VORList["vor2_MVIE_c"] = this->vor2_MVIE_c; 
  VORList["vor2_MVIE_E"] = this->vor2_MVIE_E; 
  VORList["vor2_MVIE_logdet"] = this->vor2_MVIE_logdet; 
  VORList["vor2_MVIE_volratio"] = this->vor2_MVIE_volratio; 
  VORList["vor2_MVIE_status"] = this->vor2_MVIE_status; 
  
  VORList["vor1_Dikin_E"] = this->vor1_Dikin_E; 
  VORList["vor2_Dikin_E"] = this->vor2_Dikin_E; 
  
  VORList["vis_par"] = this->vis_par; 
  
  // Flags
  VORList["isinit"] = this->isinit; 
  VORList["isset_GADJ"] = this->isset_GADJ; 
  VORList["isset_DADJ"] = this->isset_DADJ; 
  
  VORList["isset_vor2_centers"] = this->isset_vor2_centers; 
  
  VORList["isset_vor1_MVIE"] = this->isset_vor1_MVIE;
  VORList["isset_vor2_MVIE"] = this->isset_vor2_MVIE; 
  
  VORList["isset_vor1_Dikin"] = this->isset_vor1_Dikin; 
  VORList["isset_vor2_Dikin"] = this->isset_vor2_Dikin; 
  
  
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function saveRDS = base["saveRDS"];
  saveRDS(Rcpp::wrap(VORList), Rcpp::Named("file", rdsfile));
  
  return; 
}

inline void VOR::load(std::string rdsfile) {
  
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function readRDS = base["readRDS"];
  Rcpp::List invor = readRDS(Rcpp::Named("file", rdsfile));
  
  VOR_params_container par; 
  par.parallel = Rcpp::as<bool>(invor["params_parallel"]);
  par.MVIE_maxiter = Rcpp::as<int>(invor["params_MVIE_maxiter"]);
  par.MVIE_tol = Rcpp::as<double>(invor["params_MVIE_tol"]);
  par.MVIE_fix_c0 = Rcpp::as<bool>(invor["params_MVIE_fix_c0"]);
  par.seed_DADJ = Rcpp::as<bool>(invor["params_seed_DADJ"]);
  this->params = par; 
  
  // Prototypes & bounds 
  this->W = Rcpp::as<arma::mat>(invor["W"]); 
  this->nW = Rcpp::as<unsigned int>(invor["nW"]); 
  this->d = Rcpp::as<unsigned int>(invor["d"]); 
  this->lb = Rcpp::as<arma::vec>(invor["lb"]); 
  this->ub = Rcpp::as<arma::vec>(invor["ub"]); 
  this->CADJ = Rcpp::as<arma::umat>(invor["CADJ"]); 
  this->GADJ = Rcpp::as<arma::umat>(invor["GADJ"]); 
  
  this->vor1_active = Rcpp::as<arma::uvec>(invor["vor1_active"]); 
  this->vor2_active = Rcpp::as<arma::umat>(invor["vor2_active"]); 
  
  this->DADJ = Rcpp::as<arma::umat>(invor["DADJ"]); 
  
  this->vor1_centers = Rcpp::as<arma::mat>(invor["vor1_centers"]); 
  this->vor2_centers = Rcpp::as<arma::mat>(invor["vor2_centers"]);
  
  this->vor1_MVIE_c = Rcpp::as<arma::mat>(invor["vor1_MVIE_c"]); 
  this->vor1_MVIE_E = Rcpp::as<arma::cube>(invor["vor1_MVIE_E"]);
  this->vor1_MVIE_logdet = Rcpp::as<arma::vec>(invor["vor1_MVIE_logdet"]);
  this->vor1_MVIE_volratio = Rcpp::as<arma::vec>(invor["vor1_MVIE_volratio"]); 
  this->vor1_MVIE_status = Rcpp::as<arma::ivec>(invor["vor1_MVIE_status"]);
  
  this->vor2_MVIE_c = Rcpp::as<arma::mat>(invor["vor2_MVIE_c"]); 
  this->vor2_MVIE_E = Rcpp::as<arma::cube>(invor["vor2_MVIE_E"]);
  this->vor2_MVIE_logdet = Rcpp::as<arma::vec>(invor["vor2_MVIE_logdet"]);
  this->vor2_MVIE_volratio = Rcpp::as<arma::vec>(invor["vor2_MVIE_volratio"]); 
  this->vor2_MVIE_status = Rcpp::as<arma::ivec>(invor["vor2_MVIE_status"]);
  
  this->vor1_Dikin_E = Rcpp::as<arma::cube>(invor["vor1_Dikin_E"]); 
  this->vor2_Dikin_E = Rcpp::as<arma::cube>(invor["vor2_Dikin_E"]); 
  
  this->vis_par = Rcpp::as<Rcpp::List>(invor["vis_par"]);
  
  this->isinit = Rcpp::as<bool>(invor["isinit"]);
  this->isset_GADJ = Rcpp::as<bool>(invor["isset_GADJ"]);
  this->isset_DADJ = Rcpp::as<bool>(invor["isset_DADJ"]);
  this->isset_vor2_centers = Rcpp::as<bool>(invor["isset_vor2_centers"]);
  this->isset_vor1_MVIE = Rcpp::as<bool>(invor["isset_vor1_MVIE"]);
  this->isset_vor2_MVIE = Rcpp::as<bool>(invor["isset_vor2_MVIE"]);
  this->isset_vor1_Dikin = Rcpp::as<bool>(invor["isset_vor1_Dikin"]);
  this->isset_vor2_Dikin = Rcpp::as<bool>(invor["isset_vor2_Dikin"]);
  
  return; 
}
  
#endif


/* 
inline void VOR::calc_DADJ() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  
  // Perform a pre-solve
  Rcpp::Rcout << "Performing Delaunay presolve ... "; 
  VOR_DELAUNAY::Delaunay_presolve_worker presolvewkr(W); 
  if(params.parallel) presolvewkr.calc_parallel(); else presolvewkr.calc_serial(); 
  arma::umat DADJ_presolve = presolvewkr.build_DADJ(); 
  // Add CONN to the presolve list 
  DADJ_presolve += this->CADJ; 
  DADJ_presolve += this->CADJ.t(); 
  // Get list of matrix subscripts to skip 
  arma::umat skipthese = arma::ind2sub(arma::size(DADJ_presolve), arma::find(DADJ_presolve)).t(); 	
  Rcpp::Rcout << "done" << std::endl; 
  
  // Initialize a Gurobi environment
  // Sink its output to dev/null to avoid Gurobi printing to screen during parallel calls
  FILE* mystdout = stdout;
  stdout = fopen("/dev/null", "w");
  std::fclose(stdout);
  GRBEnv* env = 0;
  env = new GRBEnv();
  stdout = mystdout;
  
  // Initialize parallel structure & call it 
  VOR_DELAUNAY::Delaunay_ADJ_prlwkr wkr(env, W, lb, ub);
  wkr.set_skiplist(skipthese); 
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Clean up Gurobi environment
  delete env;
  
  this->DADJ = wkr.DADJ; 
  // add back in the skipped elements 
  this->DADJ.elem(arma::find(DADJ_presolve)).ones(); 
  this->isset_DADJ = true; 
}
*/

/*
inline arma::umat VOR::get_polydef_ADJ() {
  
  arma::umat out;
  
  if(params.polydef_ADJ == "cadj") {
    Rcpp::Rcout << "Using CADJ to build polytope definition." << std::endl;
    out = this->CADJ; 
    out += this->CADJ.t(); 
    out.elem(arma::find(out)).ones(); 
    return out; 
  }
  
  if(params.polydef_ADJ == "dadj") {
    Rcpp::Rcout << "Using DADJ to build polytope definition." << std::endl;
    if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ if parameter polydef_ADJ = 'dadj'"); 
    out = this->DADJ;
    return out; 
  } 
  
  if(params.polydef_ADJ == "none") {
    
    Rcpp::Rcout << "Unrestricted polytope definition." << std::endl;
    out = arma::ones<arma::umat>(this->nW, this->nW);
    out.diag().zeros(); 
    
    return out; 
  } 
  
  // If none of the above, error 
  Rcpp::stop("Unknown polydef_ADJ.");
  return out; 
}
 */

/*
inline arma::mat VOR::get_vor1_centers() {
  if(params.vor1_centers=="W") {
    return this->W.rows(vor1_active-1); 
  }
  
  if(params.vor1_centers=="chebyshev") {
    if(this->vor1_chebyshev_centers.n_rows != this->vor1_active.n_elem) 
      Rcpp::stop("Parameter vor1_centers='chebyshev', call $calc_vor1_chebyshev_centers first.");
    return this->vor1_chebyshev_centers; 
  }
  
  if(params.vor1_centers=="interior") {
    if(this->vor1_interior_points.n_rows != this->vor1_active.n_elem) 
      Rcpp::stop("Parameter vor1_centers='interior', call $set_vor1_interior_points first.");
    return this->vor1_interior_points; 
  }
  
  if(params.vor1_centers=="analytic") {
    if(this->vor1_analytic_centers.n_rows != this->vor1_active.n_elem) 
      Rcpp::stop("Parameter vor1_centers='analytic', call $calc_vor1_analytic_centers first.");
    return this->vor1_analytic_centers; 
  }
  
  return arma::zeros<arma::mat>(0,0); // should never get triggered, just here to satisfy compiler 
  
}

inline arma::mat VOR::get_vor2_centers() {
  
  if(params.vor2_centers=="chebyshev") {
    if(this->vor2_chebyshev_centers.n_rows != this->vor2_active.n_rows) 
      Rcpp::stop("Parameter vor2_centers='chebyshev', call $calc_vor2_chebyshev_centers first.");
    return this->vor2_chebyshev_centers; 
  }
  
  if(params.vor2_centers=="interior") {
    if(this->vor2_interior_points.n_rows != this->vor2_active.n_rows) 
      Rcpp::stop("Parameter vor2_centers='interior', call $set_vor2_interior_points first.");
    return this->vor2_interior_points; 
  }
  
  if(params.vor2_centers=="analytic") {
    if(this->vor2_analytic_centers.n_rows != this->vor2_active.n_rows) 
      Rcpp::stop("Parameter vor2_centers='analytic', call $calc_vor2_analytic_centers first.");
    return this->vor2_analytic_centers; 
  }
  
  return arma::zeros<arma::mat>(0,0); // should never get triggered, just here to satisfy compiler 
  
}
*/

/*
inline void VOR::calc_vor1_analytic_centers() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  
  // Decode which adjacency matrix to use 
  //arma::umat polydef_ADJ = this->get_polydef_ADJ(); 
  
  // Initialize the struct
  arma::uvec i_list = this->vor1_active; 
  i_list -= 1; // -1 since we're passing to c++
  //arma::mat X0 = this->get_vor1_centers(); 
  VOR_POLYINTERIOR::vor1_analytic_center_prlwkr wkr(this->W, this->DADJ, i_list, this->lb, this->ub, this->W.rows(i_list)); 
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Otherwise store the points 
  this->vor1_analytic_centers = wkr.C; 
}

inline void VOR::calc_vor2_analytic_centers() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor2_active.n_elem == 0) Rcpp::stop("Must call $set_vor2_active first");
  if(!this->isset_DADJ) Rcpp::stop("Must call $calc/set_DADJ first.");
  if(this->vor2_chebyshev_centers.n_rows != this->vor2_active.n_rows) {
    this->calc_vor2_chebyshev_centers(); 
  } 
  
  // Decode which adjacency matrix to use 
  //arma::umat polydef_ADJ = this->get_polydef_ADJ(); 
  
  
  // Initialize the struct
  arma::umat ij_list = this->vor2_active; 
  ij_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor2_analytic_center_prlwkr wkr(this->W, this->DADJ, ij_list, this->lb, this->ub, this->vor2_chebyshev_centers);
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Store the points 
  this->vor2_analytic_centers = wkr.C; 
}

inline void VOR::set_vor1_interior_points(const arma::mat& X0) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_VOR first.");
  if(this->vor1_active.n_elem == 0) Rcpp::stop("Must call $set_vor1_active first");
  
  // Decode which adjacency matrix to use 
  arma::umat polydef_ADJ = this->get_polydef_ADJ(); 
  // arma::umat polydef_ADJ; 
  // if(this->isset_DADJ) {
  //   polydef_ADJ = this->DADJ; 
  // } else {
  //   polydef_ADJ = arma::ones<arma::umat>(this->nW, this->nW); 
  //   polydef_ADJ.diag().zeros(); 
  // }
  
  
  // Initialize the struct
  arma::uvec internal_i_list = this->vor1_active; 
  internal_i_list -= 1; // -1 since we're passing to c++
  VOR_POLYINTERIOR::vor1_check_interior_prlwkr wkr(this->W, polydef_ADJ, internal_i_list, this->lb, this->ub, X0); 
  
  // Run, according to choice
  if(params.parallel) wkr.calc_parallel_progress(); else wkr.calc_serial(); 
  
  // Check if any points are not interior, if so print their indices to console and stop. 
  if(arma::any(wkr.is_interior == 0)) {
    arma::uvec not_interior = arma::find(wkr.is_interior==0); 
    not_interior += 1; 
    
    Rcpp::Rcout << "These rows of X0 are not interior to their vor2 polytopes:" << std::endl; 
    for(unsigned int i = 0; i < not_interior.n_elem; ++i) {
      Rcpp::Rcout << not_interior[i]; 
      if((i+1) % 8 == 0 || (i+1)==not_interior.n_elem) {
        Rcpp::Rcout << std::endl;
      } else {
        Rcpp::Rcout << "\t"; 
      }
    }
    
    Rcpp::stop("All points in X0 must be interior");
  }
  
  // Otherwise store the points 
  this->vor1_interior_points = X0; 
  this->isset_vor1_interior_points = true; 
}


*/
