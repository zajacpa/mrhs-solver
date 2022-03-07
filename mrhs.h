/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac 
 *      MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _MRHS_H
#define _MRHS_H


#include "mrhs.bm.h"

/// 
/// MRHS system
/// 
typedef struct {
	int nblocks; 	// number of blocks
	_bm *pM;        // left  hand side matrix blocks
	_bm *pS;        // right hand side sets
} MRHS_system;

/// Construction and destruction

/// n, m, (l, k) - fixed
MRHS_system create_mrhs_fixed(int nrows, int nblocks, int blocksize, int rhscount);

/// n, m, (l, k) as lists
MRHS_system create_mrhs_variable(int nrows, int nblocks, int blocksizes[], int rhscounts[]);

/// free memory structures
void clear_MRHS(MRHS_system *system); 

/// Random systems

void fill_mrhs_random(MRHS_system *psystem);
void fill_mrhs_random_sparse(MRHS_system *psystem);
void fill_mrhs_random_sparse_extra(MRHS_system *psystem, int density);

/// I/O

MRHS_system read_mrhs_variable(FILE *f);
int write_mrhs_variable(FILE *f, MRHS_system system);
int print_mrhs(FILE *f, MRHS_system system);

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

#endif //_MRHS_H
