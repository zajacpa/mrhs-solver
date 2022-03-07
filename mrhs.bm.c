////////////////////////////////////////////////////////////////////////
// Block  Matrix

#include <stdlib.h>
#include "mrhs.bm.h"

////////////////////////////////////////////////////////////////////////
// BM Constructors and destructors

/// Creates a dynamic BitMatrix with nrows and ncols, 
/// alloc: array of blocks
_bm create_bm(int nrows, int ncols)
{
    _bm bm; 
    //int nblocks = GET_NUM_BLOCKS(ncols);
    
    if (nrows == 0)
    {
		bm.ncols = ncols;
		bm.nrows = 0;
		bm.rows  = NULL;
		return bm;
	}
            
    bm.ncols = ncols;
    bm.nrows = nrows;
    
    bm.rows = (_block*) calloc(nrows, sizeof(_block));
    
    return bm;
}


///free space allocated to internal _bm structures
void clear_bm(_bm* pbm)
{  
	if (pbm->rows != NULL) { free(pbm->rows); }
	
    pbm->rows  = NULL;
    pbm->nrows = 0; 
}


////////////////////////////////////////////////////////////////////////////////
// Utility functions - random generation

_block random_block()
{
	return (((_block)rand()) << 32) ^ ((_block)rand() & (_block)0xffffffff);
}

///fill in with random values
void random_bm(_bm *pbm)
{
	int row;
	_block mask = BLOCK_MASK(pbm->ncols);

	for (row = 0; row < pbbm->nrows; row++)
	{     
		pbbm->rows[row] = random_block() & mask;
	}
}

///fill in with unique random values
/// PRE: 
void random_unique_bm(_bm *pbm)
{
	int filled, row;
	_block mask = BLOCK_MASK(pbm->ncols), value, tmp;
	
	for (filled = 0; filled < pbm->nrows; filled++)
	{
		//generate value, and insert sort it
		value = random_block() & mask;
		for (row = 0; row < filled; row++)
		{
			//should value be inserted here?
			if (pbm->rows[row] > value)
			{
				//if yes, shift others up
				tmp            = pbm->rows[row];
				pbm->rows[row] = value;
				value          = tmp;
			}
			//duplicate?
			else if (pbm->rows[row] == value)
			{
				//add one to value, continue with insert
				value++;
			}
		}
		pbm->rows[filled] = value;
	}
}

///fill in with random values, 
///    single    one to each column, linearly independent
void random_sparse_cols_bm(_bm *pbm)
{
	int row, col;
	_block mask = BLOCK_MASK(pbm->ncols);
     
	for (col = 0; col < pbm->ncols; col++)
	{
		mask = (ONE<<col);
		
		row = rand() % pbm->nrows;
		if (pbm->rows[row] != ZERO)
		{
			//already set to another one
			col--; 
			continue;
		}
		
		//set active variable
		pbm->rows[row] = mask;
	}   
}

/// --------------------------------------------------------------------
/// I/O
void print_block(FILE* f, _bm bm, int row)
{
	_block data = bm.rows[row]; 
	for (int bit = 0; bit < system.pM[block].ncols; bit++, data>>=1)
	{
		fprintf(f, "%llu", data&ONE);
	}
}
int read_block(FILE* f, _bm bm, int row)
{
	int c = 0;
	_block data = ZERO;
	for (bit = 0; bit < bm.ncols; bit++)
	{
		c += fscanf(f, "%1i", &read); 
		data ^= (((_block)read)<<bit); 
	}
	bm.rows[row] = data;
	return c;
}
