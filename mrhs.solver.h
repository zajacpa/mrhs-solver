/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _SOLVER_H
#define _SOLVER_H

/***************************************************************************
 * Data structures
 *
 * BlockBitMatrix : 
 *  NBLOCKS x {BLOCKSIZE}
 *  NROWS x {NBLOCKS x {BLOCKDATA}}
 *
 * BlockDiagonalBitMatrix : 
 *  NBLOCKS x {BLOCKSIZE}
 *  NBLOCKS x {NROWS x {BLOCKDATA}}
 ***************************************************************************/



//rounded up
#define GET_BL(X) (((X)+MAXBLOCKSIZE-1)/MAXBLOCKSIZE)

///main data structure for block bit matrix 
typedef struct {
   int nblocks;      // number of blocks - size of internal array in each row
   int *blocksizes;  // blocksizes in each block - max MAXBLOCKSIZE   
   int *pivots;      // number of pivots in each block (or special use)
   int nrows;        // number of rows
   int ncols;		   // total number of columns = sum of blocksizes
   _block **rows;      // storage: array of nrows pointers to arrays of blocks
} _bbm;

/// Creates a dynamic BlockBitMatrix with nrows and nblocks, 
///   each with the same blocksize, blocks filled with zeroes
/// alloc: array of pointers
/// PRE: nblocks > 0, 0 < blocksize <= MAXBLOCKSIZE, nrows > 0
_bbm* create_bbm(int nrows, int nblocks, int blocksize);

_bbm* create_bbm_new(int nrows, int nblocks, int blocksizes[]);


///free space allocated to internal _bbm structures
void free_bbm(_bbm* pbbm);

/// Creates an echelon form of matrix and computes number of pivots
///  swaps columns so that pivots form identity matrix at the start of block
/// if rhs pointer is non-zero, swaps columns also in rhs(must have same blocks)
/// NOTE: pivots are moved to MSB part, so that LSB part can be used as an index
int echelonize(_bbm *pbbm, _bbm *prhs[], _bbm **pA);


////////////////////////////////////////////////////////////////////////////////
// Tables for computation

typedef struct _te {
    _block   value;
    _block  *sm_row;
    int  first;       //first non-zero index
    struct _te *next;   
} TableEntry;

typedef struct {
    _block  mask; 
    TableEntry* *LUT;   //array of pointers
    _block* u;   
    _block  val;
    TableEntry *next;
} ActiveListEntry;

//PRE: pbbm and prhs prepared by echelonize
//TODO?: variable block sizes - this should already work 
//WORKAROUND: allows variable number of rhs by removing duplicate entries
ActiveListEntry* prepare(_bbm *pbbm, _bbm *prhs[]);

///free memory allocated to lookup tables
void free_ales(ActiveListEntry* ale, int count);

///Solver core function
typedef void (*sol_rep_fn_t)(long long int counter, _bbm *pbbm, ActiveListEntry* ale);

//front end to non-recursive call
//TODO: for multiprocessing, fork can be used and new process created for each rhs
//TODO: for threading, sol must be created for each rhs/thread 
long long int solve(ActiveListEntry* ale, _bbm *pbbm, long long int *pCount, long long int *pXors, sol_rep_fn_t report_solution);


///formula from article Ntotal
/// sum ( prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_expected(_bbm *pbbm, _bbm *prhs[]);

///formula from article Nxor
/// sum ( (m-i+1) prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_xor1(_bbm *pbbm, _bbm *prhs[]);

///formula from article Nxored
/// sum ( (1-2^(-p{i-1})) (m-i+1) prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_xor2(_bbm *pbbm, _bbm *prhs[]);

#endif //_SOLVER_H
