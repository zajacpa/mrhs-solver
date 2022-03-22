/**********************************
 * MRHS based solver
 * (C) 2016 Pavol Zajac
 *
 * I/O + testing module
 *
 * v1.3: refactoring + solution reporting + variable sized blocks
 * v1.4: minor fixes (reporting)
 * v1.6: fix input type - int to 64 bit when constructing matrix
 * v1.7: refactoring + integration of HC
 *
 * Compilation: make...
 **********************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "mrhs.bv.h"
#include "mrhs.h"
#include "mrhs.hillc.h"
#include "mrhs.solver.h"


_bbm *GlobalA = NULL;
_bv  *GlobalResults = NULL;

//TODO: create function in solver to get solution y, and to multiply y*A
int report_solution_extract_y(long long int counter, _bbm *pbbm, ActiveListEntry* ale)
{
     int block, pivot, i;
     _block value, y[pbbm->nrows];

#if (_VERBOSITY > 1)
	 fprintf(stdout, "Found solution %lli: ", counter);
#endif

     GlobalResults = (_bv*) realloc(GlobalResults, counter*(sizeof(_bv)));
     GlobalResults[counter-1] = create_bv(GlobalA->nblocks);

     //TODO: create better functions for this...
     i = 0;
     for (block = 0; block < pbbm->nblocks; block++)
     {
         value = ale[block].val;
         value>>=(pbbm->blocksizes[block]-pbbm->pivots[block]);
         for (pivot = 0; pivot < pbbm->pivots[block]; pivot++, value>>=1)
         {
             y[i++] = value&ONE;
         }
     }

     for (block = 0; block < GlobalA->nblocks; block++)
     {
         value = 0;
         for (i = 0; i < GlobalA->nrows; i++)
         {
            value ^= y[i] & GlobalA->rows[i][block];
         }
         set_bit_bv(&GlobalResults[counter-1], block, value&ONE);
#if (_VERBOSITY > 1)
		fprintf(stdout, "%01x", (unsigned) (value&ONE));
#endif
     }
#if (_VERBOSITY > 1)
	 fprintf(stdout, "\n");
#endif
	return 0;   //return 1; to find multiple solutions
}


//front end to non-recursive call
//TODO: connect with MRHS RZ solver, refactor...
long long int solve_rz(MRHS_system *system, _bv **pResults, int maxt, long long int* pTotal, long long int* pXors)
{
     ActiveListEntry* pActiveList;
     long long int count = 0;

    _bbm *pbbm, **prhs, *pA = NULL;

	//TODO: pbbm and prhs from system...
	int blocksizes[system->nblocks];
	for (int block = 0; block < system->nblocks; block++)
		blocksizes[block] = system->pM[block].ncols;

    pbbm = create_bbm_new(system->pM[0].nrows, system->nblocks, blocksizes);
    for (int block = 0; block < system->nblocks; block++)
    {
		for (int row = 0; row < pbbm->nrows; row++)
		{
			pbbm->rows[row][block] = system->pM[block].rows[row];
		}
	}

    prhs = (_bbm**) calloc(pbbm->nblocks, sizeof(_bbm*));
    for (int block = 0; block < pbbm->nblocks; block++)
    {
        prhs[block] = create_bbm(system->pS[block].nrows, 1, system->pS[0].ncols);
        for (int row = 0; row < prhs[block]->nrows; row++)
        {
			prhs[block]->rows[row][0] = system->pS[block].rows[row];
		}
	}

#if (_VERBOSITY > 1)
    int rank =
#endif
                echelonize(pbbm, prhs, &pA);

    GlobalA = pA;
#if (_VERBOSITY > 1)
	fprintf(stdout, "Starting RZ solver, system rank = %i\n", rank);
#endif

    pActiveList = prepare(pbbm, prhs);

    *pTotal = solve(pActiveList, pbbm, &count, pXors, report_solution_extract_y);

    free_ales(pActiveList, pbbm->nblocks);

    if (pA != NULL)
        free_bbm(pA);
    for (int block = 0; block < pbbm->nblocks; block++)
        free_bbm(prhs[block]);
    free(prhs);
    free_bbm(pbbm);

    *pResults = GlobalResults;

#if (_VERBOSITY > 1)
	fprintf(stdout, "RZ done\n");
#endif

   	return count;
}
