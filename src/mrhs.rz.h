/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _SOLVER_RZ_H
#define _SOLVER_RZ_H

//version 1.7

#include "mrhs.bm.h"
#include "mrhs.h"

//front end to non-recursive call
//TODO: connect with MRHS RZ solver, refactor...
long long int solve_rz(MRHS_system *system, _bv **pResults, int maxt, long long int* pCount, long long int* pXors);

#endif //_SOLVER_H
