///////////////////////////////////////////////////////////////////////
// Block Matrix

#ifndef _MRHS_BM_H
#define _MRHS_BM_H

#include <stdint.h>
#include <stdio.h>

/**********************************************************************
 * Data structure
 *
 * BitMatrix : 
 *  NROWS   x {NCOLS x {BLOCKDATA}}
 **********************************************************************/

///////////////////////////////////////////////////////////////////////
// Data Structure BM - BitMatrix

typedef uint64_t _block;

#define MAXBLOCKSIZE    64
#define ONE  ((_block)1llu)
#define ZERO ((_block)0llu)


//mask for last block
#define BLOCK_MASK(X)  ((((X) == MAXBLOCKSIZE) ? ZERO : (ONE<<(X)))-1) 

///main data structure for block bit matrix 
typedef struct {
   int nrows;          // number of rows
   int ncols;		   // number of columns
   _block *rows;      // storage: array of blocks
} _bm;

/// --------------------------------------------------------------------
/// Alloc/dealloc

/// Creates a dynamic BlockBitMatrix with nrows and ncols, 
/// alloc: array of pointers
_bm create_bm(int nrows, int ncols);

///free space allocated to internal _bbm structures
void clear_bm(_bm* pbm);

/// --------------------------------------------------------------------
/// Random data

/// random block matrix
void random_bm(_bm *pbm);

///fill in with unique random values
/// PRE: pbm->nrows < ((ONE) << pbm->ncols)
void random_unique_bm(_bm *pbm);

///fill in with random values, 
///    single    one to each column, linearly independent
void random_sparse_cols_bm(_bm *pbm);

/// --------------------------------------------------------------------
/// I/O

/// print data to f from row, sequence of 0/1
void print_block(FILE* f, _bm bm, int row);

/// read data from f from row, sequence of 0/1
int read_block(FILE* f, _bm bm, int row);

#endif //_MRHS_BBM_H
