#ifndef VORVQ_UTILS_PARALLEL_HPP
#define VORVQ_UTILS_PARALLEL_HPP

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]


/* The following re-def of "FALSE" is to fix this error. Think it has somethign to do with RcppParallel header 
 https://github.com/RcppCore/Rcpp/issues/846

RcppExports.cpp:159:5: error: no matching function for call to 'R_useDynamicSymbols'
R_useDynamicSymbols(dll, FALSE);
^~~~~~~~~~~~~~~~~~~
  /Library/Frameworks/R.framework/Resources/include/R_ext/Rdynload.h:84:10: note: candidate function not viable: no known conversion from 'int' to 'Rboolean' for 2nd argument
  Rboolean R_useDynamicSymbols(DllInfo *info, Rboolean value);
^
  1 error generated.
 */
#ifdef FALSE
#undef FALSE
#endif

#include <string>
#include <chrono> 


namespace VOR_PARALLEL {

inline unsigned int get_NumThreads_RcppParallel() {
  // Extract the environment variables from R as a list 
  Rcpp::Function Getenv("Sys.getenv");
  Rcpp::List out = Getenv(); 
  
  // If numthreads has been set via call to RcppParallel::setThreadOptions
  // the # of threads will be saved as an element of this list. 
  if(out.containsElementNamed("RCPP_PARALLEL_NUM_THREADS")) {
    std::string tmp = out["RCPP_PARALLEL_NUM_THREADS"];
    return std::stoul(tmp);
  // Otherwise, just use the default num threads - 1
  } else {
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("RcppParallel");
    Rcpp::Function Threads = pkg["defaultNumThreads"];
    return Rcpp::as<unsigned int>(Threads())-1;
  }
}

}

#endif