///////////////////////////////////////////////////////////////////////
// Block Matrix

#ifndef _MRHS_BV_H
#define _MRHS_BV_H

#include <stdint.h>
#include <stdio.h>

#include "mrhs.bm.h"

/**********************************************************************
 * Data structure
 *
 * BitVector : 
 *  sequence of bits
 **********************************************************************/


///main data structure for block bit matrix 
typedef struct {
   int ncols;		   // number of columns
   int nblocks;        // 
   _block *row;        // storage: array of blocks
} _bv;

/// --------------------------------------------------------------------
/// Alloc/dealloc

/// Creates a dynamic BlockBitMatrix with nrows and ncols, 
/// alloc: array of pointers
_bv create_bv(int ncols);

///free space allocated to internal _bbm structures
void clear_bv(_bv* pbv);

_block get_bit_bv(_bv *bv, int col);
void set_one_bv(_bv *bv, int col);
void set_zero_bv(_bv *bv, int col);
void set_bit_bv(_bv *bv, int col, _block bit);

_block proj_bv(_bv *bv, int from, int size);
void set_block_bv(_bv *bv, int from, int size, _block bit);

///swaps cols in a block
void swap_cols_bv(_bv *pv, int col1, int col2);

void print_bv(_bv *pv, FILE* f);

#endif //_MRHS_BBM_H
