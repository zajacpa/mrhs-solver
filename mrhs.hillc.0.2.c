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

////////////////////////////////////////////////////////////////////////////////
// BBM Constructors and destructors

/// Creates a dynamic BlockBitMatrix with nrows and nblocks, 
///   each with the same blocksize, blocks filled with zeroes
/// alloc: array of pointers
/// PRE: nblocks > 0, 0 < blocksize <= MAXBLOCKSIZE, nrows > 0
_bbm* create_bbm(int nrows, int nblocks, int blocksize)
{
    _bbm *pbbm; 
    int i;
    
    pbbm = (_bbm*) malloc(1 * sizeof(_bbm));
    
    pbbm->nblocks    = nblocks;
    pbbm->blocksizes = (int*) malloc(nblocks * sizeof(int));
    pbbm->pivots     = (int*) malloc(nblocks * sizeof(int));
    
    pbbm->ncols = 0;
    for (i = 0; i < nblocks; i++) {
        pbbm->blocksizes[i] = blocksize;
        pbbm->pivots[i]     = 0;
        pbbm->ncols        += blocksize;
    }

    pbbm->nrows      = nrows;
    pbbm->rows = (_block**) malloc(nrows * sizeof(_block*));
    for (i = 0; i < nrows; i++) {
        pbbm->rows[i] = (_block*) calloc(nblocks, sizeof(_block));
    }
    
    return pbbm;
}

/// Creates a dynamic BlockBitMatrix with nrows and nblocks, 
///   each with the same blocksize, blocks filled with zeroes
/// alloc: array of pointers
/// PRE: nblocks > 0, 0 < blocksize <= MAXBLOCKSIZE, nrows > 0
_bbm* create_bbm_new(int nrows, int nblocks, int blocksizes[])
{
    _bbm *pbbm; 
    int i;
    
    pbbm = (_bbm*) malloc(1 * sizeof(_bbm));
    
    pbbm->nblocks    = nblocks;
    pbbm->blocksizes = (int*) malloc(nblocks * sizeof(int));
    pbbm->pivots     = (int*) malloc(nblocks * sizeof(int));
    
    pbbm->ncols = 0;
    for (i = 0; i < nblocks; i++) {
        pbbm->blocksizes[i] = blocksizes[i];
        pbbm->pivots[i]     = 0;
        pbbm->ncols        += blocksizes[i];
    }

    pbbm->nrows      = nrows;
    pbbm->rows = (_block**) malloc(nrows * sizeof(_block*));
    for (i = 0; i < nrows; i++) {
        pbbm->rows[i] = (_block*) calloc(nblocks, sizeof(_block));
    }
    
    return pbbm;
}

///free space allocated to internal _bbm structures
void free_bbm(_bbm* pbbm)
{
    int i; 
    for (i = 0; i < pbbm->nrows; i++) {
        free(pbbm->rows[i]);
    }    
    free(pbbm->rows);
   
    free(pbbm->pivots);
    free(pbbm->blocksizes);

    free(pbbm);
}

////////////////////////////////////////////////////////////////////////////////
// Sparse bit array

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

SparseBitArray* to_sparse_bit_array(_bbm *pbbm)
{
	//ASSERT(pbbm->nblocks == 1)
	if (pbbm->nblocks < 1 || pbbm->nrows < 1)
		return NULL;
	
	SparseBitArray* sba = (SparseBitArray*) calloc(1, sizeof(SparseBitArray));
	int block = 0;
	int width = pbbm->blocksizes[block];
	
	//simple algorithm: only one value, place it at given offset
	if (pbbm->nrows == 1)
	{
		sba->offset = pbbm->rows[0][block];
		sba->value  = ONE;
		sba->next   = NULL;
		return sba;
	}
	//simple algorithm 2, all values within one block
	if ((ONE << width) <= MAXBLOCKSIZE)
	{
		sba->offset = ZERO;
		sba->value  = ZERO;
		
		for (int row = 0; row < pbbm->nrows; row++)
		{
			sba->value ^= (ONE << pbbm->rows[row][block]);
		}
		
		sba->next = NULL;
		return sba;
	}
	
	//complex algorithm, multiple bit blocks
	
	//step1: sort values
	_block values[pbbm->nrows];
	for (int row = 0; row < pbbm->nrows; row++)
	{
		values[row] = pbbm->rows[row][block];
	}
	qsort(values, pbbm->nrows, sizeof(_block), cmp_block);
	
	//step2: initial position
	SparseBitArray* start = sba;
	sba->offset = values[0];
	sba->value  = ONE;
	sba->next   = NULL;
		
	for (int row = 1; row < pbbm->nrows; row++)
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
   _bbm *pbbm;
   SparseBitArray **rhs; // sparse bit array for each RHS set
} CompressedMRHS;


void add_row(_block out[], int row, CompressedMRHS* cmrhs)
{
	for (int block = 0; block < cmrhs->pbbm->nblocks; block++)
	{
		out[block] ^= cmrhs->pbbm->rows[row][block];
	}		
}

//PRE: pbbm and prhs are valid MRHS system
CompressedMRHS* prepare(_bbm *pbbm, _bbm *prhs[])
{   
    CompressedMRHS *cmrhs;
    
    //allocate compressed representation of MRHS
    cmrhs = (CompressedMRHS*) calloc(1, sizeof(CompressedMRHS));
    cmrhs->pbbm = pbbm;
      
    cmrhs->rhs = (SparseBitArray**) calloc(cmrhs->pbbm->nblocks, sizeof(SparseBitArray*));
    for (int block = 0; block < cmrhs->pbbm->nblocks; block++)
    {
		cmrhs->rhs[block] = to_sparse_bit_array(prhs[block]);
	}
    
    return cmrhs;
}


void free_cmrhs(CompressedMRHS* cmrhs)
{
    int row, block;
    
    for (block = 0; block < cmrhs->pbbm->nblocks; block++)
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
    for (int block = 0; block < cmrhs->pbbm->nblocks; block++)
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
long long int solve(_bbm *pbbm, _bbm *prhs[], int maxs, long long int* pCount, long long int* pRestarts)
{
	CompressedMRHS* cmrhs = prepare(pbbm, prhs);
	long long int count = 0;
	int diff, bestdiff, bestix, restart;
	
	_block *solution = (_block*) malloc(pbbm->nrows * sizeof(_block));
	
	_block *rhs  =(_block*) malloc(cmrhs->pbbm->nblocks * sizeof(_block));
	_block *tmp  =(_block*) malloc(cmrhs->pbbm->nblocks * sizeof(_block));
	
	for (restart = 0; clock() < maxs; restart++)
	{	
		//initialize random solution
		memset(rhs, 0, cmrhs->pbbm->nblocks * sizeof(_block));
		for (int row = 0; row < cmrhs->pbbm->nrows; row++)
		{
			solution[row] = rand() % 2;
			if (solution[row] != 0)
			{
				add_row(rhs, row, cmrhs);
			}
		}
		
		//check solution
		bestdiff = evaluate(rhs, cmrhs);
		while (bestdiff > 0)
		{
			bestix = -1;
			//for each row:
			for (int row = 0; row < pbbm->nrows; row++)
			{
				memcpy(tmp, rhs, cmrhs->pbbm->nblocks * sizeof(_block));
				add_row(tmp, row, cmrhs);
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
				add_row(rhs, bestix, cmrhs);
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
		printf("Solution found in %i restarts: ", restart);
		for (int row = 0; row < cmrhs->pbbm->nrows; row++)
		{
			printf("%lu", solution[row]);
		}
		printf("\n");
	#endif    
		
		return 1;
	}
	else
	{
		return 0;
	}
}

