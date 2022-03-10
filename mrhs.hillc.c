/**********************************
 * MRHS based solver
 * (C) 2016 Pavol Zajac  
 * 
 * library file
 *
 * v0.1: initial algorithm
 **********************************/

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>

#include "mrhs.hillc.h"

/// ////////////////////////////////////////////////////////////////////
/// Sparse bit array

typedef struct _sba {
	_block offset;
	_block value;
	struct _sba *next;
} SparseBitArray;

//get value in sparse bit array (ZERO or ONE)
_block sba_value_at(SparseBitArray* sba, _block position)
{
	while (sba != NULL && sba->offset + MAXBLOCKSIZE <= position)
		sba = sba->next;
	
	//after the end, all zeroes
	if (sba == NULL || sba->offset > position)
		return ZERO;
	
	//in the block...
	return ((sba->value)>>(position - sba->offset)) & ONE;
} 

int cmp_block(const void* a, const void* b)
{
	if (*(_block*)a == *(_block*)b) return 0;
	if (*(_block*)a >  *(_block*)b) return 1;
	return -1;
}

SparseBitArray* to_sparse_bit_array(_bm bm)
{
	//ASSERT(pbbm->nblocks == 1)
	if (bm.nrows < 1)
		return NULL;
	
	SparseBitArray* sba = (SparseBitArray*) calloc(1, sizeof(SparseBitArray));
	int width = bm.ncols;
	
	//simple algorithm: only one value, place it at given offset
	if (bm.nrows == 1)
	{
		sba->offset = bm.rows[0];
		sba->value  = ONE;
		sba->next   = NULL;
		return sba;
	}
	//simple algorithm 2, all values within one block
	if ((ONE << width) <= MAXBLOCKSIZE)
	{
		sba->offset = ZERO;
		sba->value  = ZERO;
		
		for (int row = 0; row < bm.nrows; row++)
		{
			sba->value ^= (ONE << bm.rows[row]);
		}
		
		sba->next = NULL;
		return sba;
	}
	
	//complex algorithm, multiple bit blocks
	
	//step1: sort values
	_block values[bm.nrows];
	for (int row = 0; row < bm.nrows; row++)
	{
		values[row] = bm.rows[row];
	}
	qsort(values, bm.nrows, sizeof(_block), cmp_block);
	
	//step2: initial position
	SparseBitArray* start = sba;
	sba->offset = values[0];
	sba->value  = ONE;
	sba->next   = NULL;
		
	for (int row = 1; row < bm.nrows; row++)
	{
		//add other positions
		// if within continues block, store it
		if (values[row] - sba->offset < MAXBLOCKSIZE)
		{
			sba->value ^= (ONE << (values[row] - sba->offset));
		}
		//otherwise start new block
		else
		{
			sba->next = (SparseBitArray*) calloc(1, sizeof(SparseBitArray));
			
			sba = sba->next;
			sba->next = NULL;
			sba->offset = values[row];
			sba->value = ONE;
		}
	}
		
	return start;	
}


void free_sba(SparseBitArray* sba)
{
	while (sba != NULL)
	{
		SparseBitArray* tmp = sba->next;
		free(sba);
		sba = tmp;
	}
}

////////////////////////////////////////////////////////////////////////////////
// MRHS representation:
//   rows     as bit arrays
//   RHS sets as sparse bit arrays

typedef struct _cmrhs {
   int  nblocks;
   const _bm* pM;
   SparseBitArray **rhs; // sparse bit array for each RHS set
} CompressedMRHS;


void add_row_hc(_block out[], int row, CompressedMRHS* cmrhs)
{
	for (int block = 0; block < cmrhs->nblocks; block++)
	{
		out[block] ^= cmrhs->pM[block].rows[row];
	}		
}

//PRE: pbbm and prhs are valid MRHS system
CompressedMRHS* prepare_hc(MRHS_system system)
{   
    CompressedMRHS *cmrhs;
    
    //allocate compressed representation of MRHS
    cmrhs = (CompressedMRHS*) calloc(1, sizeof(CompressedMRHS));
    cmrhs->nblocks = system.nblocks;
    cmrhs->pM      = system.pM;
    
    cmrhs->rhs = (SparseBitArray**) calloc(cmrhs->nblocks, sizeof(SparseBitArray*));
    for (int block = 0; block < cmrhs->nblocks; block++)
    {
		cmrhs->rhs[block] = to_sparse_bit_array(system.pS[block]);
	}
    
    return cmrhs;
}


void free_cmrhs(CompressedMRHS* cmrhs)
{
    int row, block;
    
    for (block = 0; block < cmrhs->nblocks; block++)
    {
		free_sba(cmrhs->rhs[block]);
	}
    free(cmrhs->rhs);
    
    free(cmrhs);	
}

int evaluate(_block rhs[], CompressedMRHS* cmrhs)
{
	_block value;
	_block sum = 0;
    for (int block = 0; block < cmrhs->nblocks; block++)
	{
		value = sba_value_at(cmrhs->rhs[block], rhs[block]);
		sum += (ONE-value);  //if not in RHS, add one to evaluate
	}
	return (int) sum;
}

#if (_VERBOSITY > 1)
#include <stdio.h>
#endif    

//PRE: pbbm and prhs are valid MRHS system
long long int solve_hc(MRHS_system system, _bv **pResults, int maxt, long long int* pCount, long long int* pRestarts)
{
	int retval;
	if (system.nblocks == 0)
		return 0;
	
	CompressedMRHS* cmrhs = prepare_hc(system);
	long long int count = 0;
	int diff, bestdiff, bestix, restart, nrows;
	
	nrows = system.pM[0].nrows;
	_block *solution = (_block*) malloc( nrows * sizeof(_block));
	
	_block *rhs  =(_block*) malloc(cmrhs->nblocks * sizeof(_block));
	_block *tmp  =(_block*) malloc(cmrhs->nblocks * sizeof(_block));
	
	for (restart = 0; clock() < maxt; restart++)
	{	
		//initialize random solution
		memset(rhs, 0, cmrhs->nblocks * sizeof(_block));
		for (int row = 0; row < nrows; row++)
		{
			solution[row] = rand() % 2;
			if (solution[row] != 0)
			{
				add_row_hc(rhs, row, cmrhs);
			}
		}
		
		//check solution
		bestdiff = evaluate(rhs, cmrhs);
		while (bestdiff > 0)
		{
			bestix = -1;
			//for each row:
			for (int row = 0; row < nrows; row++)
			{
				memcpy(tmp, rhs, cmrhs->nblocks * sizeof(_block));
				add_row_hc(tmp, row, cmrhs);
				diff = evaluate(tmp, cmrhs);
				
				//evaluate, store if better
				if (diff < bestdiff)
				{
					bestdiff = diff;
					bestix   = row;
				}
				count++;
				if (bestdiff == 0)
					break;
			}
			//change to better
			if (bestix > -1)
			{
				add_row_hc(rhs, bestix, cmrhs);
				solution[bestix] ^= ONE;
			}
			//no change to better
			else
			{
				break; //while bestdiff > 0...
			}
		}
		//solution found, do not restart
		if (bestdiff == 0)
		{
			break;
		}
	}
	//solution found?
	*pCount    = count;
	*pRestarts = restart;
	
	
	if (bestdiff == 0)
	{
	#if (_VERBOSITY > 1)
		fprintf(stdout, "Solution found in %i restarts\n", restart);
	#endif    

		*pResults  = (_bv*) malloc(sizeof(_bv));
		**pResults = create_bv(nrows);
		for (int row = 0; row < nrows; row++)
		{
			set_bit_bv(*pResults, row, solution[row]);
		}
		
		retval = 1;
	}
	else
	{
		retval = 0;
	}
	
	free_cmrhs(cmrhs);
	free(solution);
	free(rhs);
	free(tmp);
	
	return retval;
}

