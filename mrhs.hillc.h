/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _SOLVER_H
#define _SOLVER_H

#include <stdint.h>

#include "mrhs.bm.h"


long long int solve(_bbm *pbbm, _bbm *prhs[], int maxr, long long int* pCount, long long int* pRestarts);

#endif //_SOLVER_H
