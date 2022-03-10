////////////////////////////////////////////////////////////////////////
// Block  Matrix

#include <stdlib.h>
#include <stdio.h>
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


_block get_bit_bm(_bm *bm, int row, int col)
{
	return (bm->rows[row] >> col) & ONE;
}
void set_one_bm(_bm *bm, int row, int col)
{
	bm->rows[row] |= (ONE << col);	
}
void set_zero_bm(_bm *bm, int row, int col)
{
	bm->rows[row] &= ~(ONE << col);	
}
void set_bit_bm(_bm *bm, int row, int col, _block bit)
{
	bm->rows[row] &= ~(ONE << col);	
	bm->rows[row] |= ((bit&ONE) << col);	
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

	for (row = 0; row < pbm->nrows; row++)
	{     
		pbm->rows[row] = random_block() & mask;
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
	_block mask;
     
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
void print_block_bm(FILE* f, _bm bm, int row)
{
	_block data = bm.rows[row]; 
	for (int bit = 0; bit < bm.ncols; bit++, data>>=1)
	{
		fprintf(f, "%lu", data&ONE);
	}
}
int read_block_bm(FILE* f, _bm bm, int row)
{
	int c = 0, read;
	_block data = ZERO;
	for (int bit = 0; bit < bm.ncols; bit++)
	{
		c += fscanf(f, "%1i", &read); 
		data ^= (((_block)read)<<bit); 
	}
	bm.rows[row] = data;
	return c;
}

/// --------------------------------------------------------------------
/// Linear algebra

///row dest ^= source 
void add_row_bm(_bm *pbm, int dest, int source)
{
     pbm->rows[dest] ^= pbm->rows[source];
}

///swaps rows i and j
void swap_row_bm(_bm *pbm, int i, int j)
{
     _block tmp   = pbm->rows[i]; 
     pbm->rows[i] = pbm->rows[j];
     pbm->rows[j] = tmp;
}

///swaps cols in a block
void swap_cols_bm(_bm *pbm, int col1, int col2)
{
     _block data;
     _block mask = ((ONE)<<col1) | ((ONE)<<col2);  
     
     for (int row = 0; row < pbm->nrows; row++)
     {
         data  = pbm->rows[row];
         //if both are same, no change, otherwise change both of them
         data ^= ( ((data >> col1)^(data>>col2))&ONE )*mask; //constant ONE
         
         pbm->rows[row] = data;
     }
}

///find first non-zero element in given block/mask, from given row index 
/// -1 if no such exists
int find_pivot_bm(_bm *pbm, int col, int from)
{
	_block mask = (ONE)<<col;
    while (from < pbm->nrows)
    {
        if ((pbm->rows[from] & mask) == ONE)
            return from;
        from++;
    }
    return -1;
}
