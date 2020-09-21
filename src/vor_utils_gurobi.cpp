/* This file exists solely to initialize a Gurobi environment 
 * into the user's R environment. This Gurobi environment will 
 * be used any time a function is called which relies on the Gurobi
 * solver (calculating DADJ, computing Chebyshev centers)
 * 
 */
#include "vor_utils_gurobi.hpp"

const GRBEnv* GUROBIENV = new GRBEnv();  

