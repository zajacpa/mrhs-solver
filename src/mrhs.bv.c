////////////////////////////////////////////////////////////////////////
// Block  Matrix

#include <stdlib.h>
#include <stdio.h>

#include "mrhs.bm.h"
#include "mrhs.bv.h"

////////////////////////////////////////////////////////////////////////
// BV Constructors and destructors


/// Creates a dynamic BlockBitMatrix with nrows and ncols,
/// alloc: array of pointers
_bv create_bv(int ncols)
{
    _bv bv;

    if (ncols == 0)
    {
		bv.ncols = ncols;
		bv.nblocks = 0;
		bv.row  = NULL;
		return bv;
	}

    int nblocks = GET_NUM_BLOCKS(ncols);
    bv.ncols   = ncols;
    bv.nblocks = nblocks;

    bv.row = (_block*) calloc(nblocks, sizeof(_block));

    return bv;
}

///free space allocated to internal _bbm structures
void clear_bv(_bv* pbv)
{
	if (pbv->row != NULL) { free(pbv->row); }

    pbv->row     = NULL;
    pbv->nblocks = 0;
    pbv->ncols   = 0;
}

_block get_bit_bv(const _bv *bv, int col)
{
	return (bv->row[col / MAXBLOCKSIZE] >> (col % MAXBLOCKSIZE)) & ONE;
}

void set_one_bv(_bv *bv, int col)
{
	int block = col / MAXBLOCKSIZE;
	col %= MAXBLOCKSIZE;
	bv->row[block] |= (ONE << col);
}
void set_zero_bv(_bv *bv, int col)
{
	int block = col / MAXBLOCKSIZE;
	col %= MAXBLOCKSIZE;
	bv->row[block] &= ~(ONE << col);
}
void set_bit_bv(_bv *bv, int col, _block bit)
{
	int block = col / MAXBLOCKSIZE;
	col %= MAXBLOCKSIZE;
	bv->row[block] &= ~(ONE << col);
	bv->row[block] |= ((bit&ONE) << col);
}

//PRE: size < BLOCKSIZE
_block proj_bv(_bv *bv, int from, int size)
{
	int block   = from / MAXBLOCKSIZE;
	int offset  = from % MAXBLOCKSIZE;
	_block bits = (bv->row[block] >> offset);

	//fill in bits from next block
	bits ^= (bv->row[block+1] << (MAXBLOCKSIZE-offset));

	//mask out the result
	return bits & BLOCK_MASK(size);
}
void set_block_bv(_bv *bv, int from, int size, _block bits)
{
	int block1  = from / MAXBLOCKSIZE;
	int block2  = (from+size) / MAXBLOCKSIZE;
	int start   = from % MAXBLOCKSIZE;
	int end     = (from+size) % MAXBLOCKSIZE;

	if (block2 > block1)
	{
		//upper part is in next block
		//set upper part
		bv->row[block2] &= ~BLOCK_MASK(end);
		bv->row[block2] ^= (bits >> (MAXBLOCKSIZE-start))&BLOCK_MASK(end);
		//set lower part
		bv->row[block1] &= BLOCK_MASK(start);
		bv->row[block1] ^= (bits << (start));
	}
	else
	{
		//upper part is in the same block
		bits &= BLOCK_MASK(size);   //zero out possible excess bits
		//set lower part: mask out replaced bits, add new bits
		bv->row[block1] &= ~(BLOCK_MASK(start)^BLOCK_MASK(end));
		bv->row[block1] ^= (bits << (start));
	}
}

///swaps cols in a block
void swap_cols_bv(_bv *bv, int col1, int col2)
{
	_block bit1 = get_bit_bv(bv, col1);
	_block bit2 = get_bit_bv(bv, col2);
	set_bit_bv(bv, col1, bit2);
	set_bit_bv(bv, col2, bit1);
}

///find first non-zero position >= given startpos ("pivot" in linear algebra)
/// returns -1 if not found
int find_nonzero(_bv *bv, int startpos)
{
    //TODO: more efficient version for sparse bv...
    while (startpos < bv->ncols)
    {
        if (get_bit_bv(bv, startpos) == ONE)
            return startpos;
        startpos++;
    }
    return -1;
}

_block is_non_zero_bv(_bv *bv)
{
    _block sum = 0;
    for (int block = 0; block < bv->nblocks; block++)
        sum |= bv->row[block];
    return sum == ZERO ? ZERO : ONE;
}

//bv1 X= bv2
void and_bv(_bv *bv1, _bv *bv2)
{
    for (int block = 0; block < bv1->nblocks; block++)
        bv1->row[block] &= bv2->row[block];
}
void or_bv(_bv *bv1, _bv *bv2)
{
    for (int block = 0; block < bv1->nblocks; block++)
        bv1->row[block] |= bv2->row[block];
}
void xor_bv(_bv *bv1, _bv *bv2)
{
    for (int block = 0; block < bv1->nblocks; block++)
        bv1->row[block] ^= bv2->row[block];
}
void inv_bv(_bv *bv1)
{
    for (int block = 0; block < bv1->nblocks; block++)
        bv1->row[block] = ~bv1->row[block];
}


///print bit vector to file
void print_bv(_bv *bv, FILE* f)
{
	int block = 0;
	for (block = 0; block < bv->nblocks-1; block++)
	{
		for (int bit = 0; bit < MAXBLOCKSIZE; bit++)
			fprintf(f, "%01x", (unsigned int) ((bv->row[block] >> bit)&ONE));
	}
	for (int bit = 0; bit < LASTBLOCKSIZE(bv->ncols); bit++)
		fprintf(f, "%01x", (unsigned int) ((bv->row[block] >> bit)&ONE));
}

///print bit vector to file
void random_bv(_bv *bv)
{
	int block = 0;
	_block bits;
	for (block = 0; block < bv->nblocks-1; block++)
	{
		bits = (rand() & 0xffffull) ^
                         ((rand() & 0xffffull) << 16) ^
                         ((rand() & 0xffffull) << 32) ^
                         ((rand() & 0xffffull) << 48);
        bv->row[block] = bits;
	}
    bits = (rand() & 0xffffull) ^
                         ((rand() & 0xffffull) << 16) ^
                         ((rand() & 0xffffull) << 32) ^
                         ((rand() & 0xffffull) << 48);
    bv->row[block] = bits & BLOCK_MASK(LASTBLOCKSIZE(bv->ncols));
}
