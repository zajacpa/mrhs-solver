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


/// --------------------------------------------------------------------
/// Alloc/dealloc

/// Creates a dynamic BlockBitMatrix with nrows and ncols,
/// alloc: array of pointers
_bv create_bv(int ncols);

///free space allocated to internal _bbm structures
void clear_bv(_bv* pbv);

///random bit vector
void random_bv(_bv* pbv);


_block is_non_zero_bv(_bv *bv);

//bv1 X= bv2
void and_bv(_bv *bv1, _bv *bv2);
void or_bv(_bv *bv1, _bv *bv2);
void xor_bv(_bv *bv1, _bv *bv2);
void inv_bv(_bv *bv1);


_block get_bit_bv(const _bv *bv, int col);
void set_one_bv(_bv *bv, int col);
void set_zero_bv(_bv *bv, int col);
void set_bit_bv(_bv *bv, int col, _block bit);

_block proj_bv(_bv *bv, int from, int size);
void set_block_bv(_bv *bv, int from, int size, _block bit);

///swaps cols in a block
void swap_cols_bv(_bv *pv, int col1, int col2);

///find first non-zero position >= given startpos ("pivot" in linear algebra)
/// returns -1 if not found
int find_nonzero(_bv *bv, int startpos);

void print_bv(_bv *pv, FILE* f);

#endif //_MRHS_BBM_H
