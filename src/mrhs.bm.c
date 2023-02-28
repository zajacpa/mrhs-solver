////////////////////////////////////////////////////////////////////////
// Block  Matrix

#include <stdlib.h>
#include <stdio.h>
#include "mrhs.bm.h"
#include "mrhs.bv.h"

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
                //add one to value, continue with insert (cyclic: 111->000)
                if (value == mask)
                {
                    value = 0;
                    row = -1;
                }
                else
                {
                    value++;
                }
            }
        }
        pbm->rows[filled] = value;
    }
}

///fill in pbm based on AND gate + random constant
/// PRE: nrows = 4, ncols = 3
void random_and_bm(_bm *pbm)
{
    if (pbm->nrows != 4 || pbm->ncols != 3)
        return;

    _block constant = random_block() & BLOCK_MASK(pbm->ncols); 

    pbm->rows[0] = 0x0 ^ constant;  //000 + c
    pbm->rows[1] = 0x1 ^ constant;  //001 + c
    pbm->rows[2] = 0x2 ^ constant;  //010 + c
    pbm->rows[3] = 0x7 ^ constant;  //111 + c
}

///fill in with random values for AND inputs, and single one for AND output,
/// PRE: ncols = 3, output_row < nrows
void random_and_cols_bm(_bm *pbm, int output_row)
{
    if (pbm->nrows < output_row || output_row < 0 || pbm->ncols != 3)
        return;

    //random value in first 2 columns
    _block mask = BLOCK_MASK(pbm->ncols-1);
    int row;

	for (row = 0; row < output_row; row++)
	{
		pbm->rows[row] = random_block() & mask;
	}

    //special output_row
    pbm->rows[output_row] = (ONE << 2);
    
}

///fill in with random values for AND inputs, and single one for AND output,
/// PRE: ncols = 3, output_row < nrows
void random_sparse_and_cols_bm(_bm *pbm, int output_row, int density)
{
    int special = 0, row, col;
    _block mask;
    
    if (output_row < 0 || pbm->ncols != 3)
        return;
    
    if (output_row < pbm->nrows)
    {
        //special output_row / column
        pbm->rows[output_row] = (ONE << 2);
        special = 1;
    }
    else
    {
        output_row = pbm->nrows;
    }
    
	for (col = 0; col < pbm->ncols - special; col++)
	{
		mask = (ONE<<col);

        //set active variable 1
        for (int i = 0; i <= density; i++)
        {
            row = rand() % output_row;        
            pbm->rows[row] |= mask;
        }
	}    
}

///fill in with random values,
///    single    one to each column, linearly independent
/// for correct lin independence PRE: pbm->nrows >> pbm->ncols
void random_sparse_cols_bm(_bm *pbm)
{
	int row, col;
	_block mask;

	for (col = 0; col < pbm->ncols; col++)
	{
		mask = (ONE<<col);

		row = rand() % pbm->nrows;
		//after && -> failsafe for system with low number of unknowns
		if (pbm->rows[row] != ZERO && pbm->ncols < pbm->nrows)
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
int print_block_bm(FILE* f, _bm bm, int row)
{
	_block data = bm.rows[row];
	int sum = 0;
	for (int bit = 0; bit < bm.ncols; bit++, data>>=1)
	{
		sum += fprintf(f, "%lu", data&ONE);
	}
	return sum;
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

///transform column of block matrix to bit vector
_bv get_column_bm(_bm* pbm, int col)
{
    _bv bv = create_bv(pbm->nrows);
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (get_bit_bm(pbm, row, col) == ONE)
            set_one_bv(&bv, row);
    }
    return bv;
}

///get bit vector indicating active rows
_bv get_active_rows_bm(_bm* pbm)
{
    _bv bv = create_bv(pbm->nrows);
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (pbm->rows[row] != ZERO)
            set_one_bv(&bv, row);
    }
    return bv;
}

///remove rows indicated by zeroes in rowmask
int remove_rows_bm(_bm* pbm, _bv* rowmask)
{
    int insert = 0;
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (get_bit_bv(rowmask, row) == ONE)
        {
            //keep this row
            pbm->rows[insert] = pbm->rows[row];
            insert++;
        }
    }
    pbm->nrows = insert;
    return pbm->nrows;
}

///transform column of block matrix to bit vector
void add_column_bm(_bm* pbm, _bv *column, int col)
{
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (get_bit_bv(column, row) == ONE)
            pbm->rows[row] ^= (ONE<<col);
    }
}

void add_constant_bm(_bm* pbm, _block c, int col)
{
    if (c == ZERO)
        return;

    for (int row = 0; row < pbm->nrows; row++)
    {
            pbm->rows[row] ^= (ONE<<col);
    }
}

//vector times matrix equals vector:
//PRE: correct dimensions of bit vector / matrix
_block multiply_bv_x_bm(const _bv* pbv, const _bm* pbm)
{
    _block out = ZERO;
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (get_bit_bv(pbv, row) == ONE)
            out ^= pbm->rows[row];
    }
    return out;
}

//unsorted search of x in block
int index_of_block_in_bm(const _bm* pbm, _block x)
{
    x &= BLOCK_MASK(pbm->ncols);
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (x == pbm->rows[row])
            return row;
    }
    return -1;
}

//check whether x is in bm, if it is not, replace nearest (lower valued) block
_block ensure_block_in_bm(_bm* pbm, _block x)
{
    int pos = 0;
    x &= BLOCK_MASK(pbm->ncols);
    for (int row = 0; row < pbm->nrows; row++)
    {
        if (x == pbm->rows[row])
        {
            pos = row;
            break;
        }
        if (x > pbm->rows[row])
        {
            pos = row;
        }
    }
    _block out = pbm->rows[pos];
    pbm->rows[pos] = x;
    return out;
}

