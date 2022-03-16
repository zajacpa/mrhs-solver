/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _SOLVER_HC_H
#define _SOLVER_HC_H

#include <stdint.h>

#include "mrhs.bm.h"
#include "mrhs.bv.h"
#include "mrhs.h"


long long int solve_hc(MRHS_system *system, _bv **pResults, int maxt, long long int* pCount, long long int* pRestarts);

#endif //_SOLVER_H
