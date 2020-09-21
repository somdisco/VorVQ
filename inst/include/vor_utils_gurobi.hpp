#ifndef VORVQ_UTILS_GUROBI_HPP
#define VORVQ_UTILS_GUROBI_HPP

#include "gurobi_c++.h"
extern const GRBEnv* GUROBIENV;

// Initialize a Gurobi environment and prevent the "Academic License" from printing to console 
// FILE* mystdout = stdout;
// stdout = fopen("/dev/null", "w");
// std::fclose(stdout);
// GRBEnv* env = 0;
// env = new GRBEnv();
// stdout = mystdout;


#endif