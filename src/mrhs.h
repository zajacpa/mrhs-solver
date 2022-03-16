/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac
 *      MRHS Solver Based on Linear Algebra and Exhaustive Search
 */

#ifndef _MRHS_H
#define _MRHS_H

#include <stdio.h>
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

///substitute given linear equation into system
int linear_substitution(MRHS_system *system, _bv *column, _block rhs);
///remove linear equations from the system
int remove_linear(MRHS_system *system);
int remove_empty(MRHS_system *system);

#endif //_MRHS_H
