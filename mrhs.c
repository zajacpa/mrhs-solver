/***
 * MRHS solver interface
 * See: Håvard Raddum and Pavol Zajac 
 *      MRHS Solver Based on Linear Algebra and Exhaustive Search
 */
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "mrhs.bm.h"
#include "mrhs.h"


/// n, m, (l, k) - fixed
MRHS_system create_mrhs_fixed(int nrows, int nblocks, int blocksize, int rhscount)
{
	int blocksizes[nblocks];
	int rhscounts[nblocks];
	memset(blocksizes, blocksize, nblocks * sizeof(blocksizes[0]));
	memset( rhscounts,  rhscount, nblocks * sizeof( rhscounts[0]));
	
	return create_mrhs_variable(nrows, nblocks, blocksizes, rhscounts);
}

/// n, m, (l, k) as lists
MRHS_system create_mrhs_variable(int nrows, int nblocks, int blocksizes[], int rhscounts[])
{
	MRHS_system system;
	
	if (nblocks == 0)
	{
		system.nblocks = 0;
		system.pM      = NULL;
		system.pS      = NULL;
		return system;
	}
	
	//else, non-zero blocks
	system.nblocks = nblocks;
	system.pM      = calloc(nblocks, sizeof(_bm));
	system.pS      = calloc(nblocks, sizeof(_bm));
	for (int block = 0; block < psystem->nblocks; block++)
	{
		system.pM[block] = create_bm(nrows, blocksizes[block]);
		system.pS[block] = create_bm(rhscounts[block], blocksizes[block]);
	}
	return system;
}

/// release memory allocated for MRHS system
void clear_MRHS(MRHS_system *psystem)
{
	if (psystem->nblocks == 0)	//ASSERT ...
		return;
	
	for (int block = 0; block < psystem->nblocks; block++)
	{
		clear_bm(psystem->pM);
		clear_bm(psystem->pS);
	}
	free(psystem->pM);
	free(psystem->pS);

	psystem->nblocks = 0;
	psystem->pM      = NULL;
	psystem->pS      = NULL;	
} 

/// ///////////////////////////////////////////////////////////////////

/// Fill MRHS system with random data
void fill_mrhs_random(MRHS_system *psystem)
{
	for (int block = 0; block < psystem->nblocks; block++)
	{
		random_bm(psystem->pM);
		random_unique_bm(psystem->pS);
	}	
}

/// Fill MRHS system with random data, 
///   single one in each linearly independent column
void fill_mrhs_random_sparse(MRHS_system *psystem)
{
	for (int block = 0; block < psystem->nblocks; block++)
	{
		random_sparse_cols_bm(psystem->pM);
		random_unique_bm(psystem->pS);
	}	
}

/// Fill MRHS system with random data
///  M is sparse -> one 1 in each column + density number of ones
void fill_mrhs_random_sparse_extra(MRHS_system *psystem, int density)
{
	fill_mrhs_random_sparse(psystem);
    for (i = 0; i < density; i++)
    {
		block = rand() % psystem->nblocks;
		row   = rand() % psystem->pM[block].nrows;
		col   = rand() % psystem->pM[block].ncols;
		psystem->pM[block].rows[row] |= (ONE << col);     	   
	}	
}

/// I/O

/// deserialize system
MRHS_system read_mrhs_variable(FILE *f)
{
	MRHS_system system;
	char c;
	int row, bit, block, read, n, m;
	int *k, *l;

	//read header
	fscanf(f, "%i", &n);
	fscanf(f, "%i", &m);

	k = (int*)calloc(m, sizeof(int));
	l = (int*)calloc(m, sizeof(int));
	for (i = 0; i < m; i++)
	{
	 fscanf(f, "%i", &l[i]);
	 fscanf(f, "%i", &k[i]);
	}

	//prepare system
	system = create_mrhs_variable(n, m, l, k);

	//no longer needed
    free(l);
    free(k);

	//read rows of M
	for (row = 0; row < n; row++)
	{   
		// skip until beggining of row is found
		while (fscanf(f, "%c", &c) && c != '[') { /* skip */ } 
			
		for (block = 0; block < system.nblocks; block++)
		{
			read_block(f, system.pM[block], row);
		}   
		//remove terminator
		fscanf(f, "]\n"); 
	}

	//read blocks of S
	for (block = 0; block < system.nblocks; block++)
	{   
		for (row = 0; row < system.pS[block].nrows; row++)
		{
			// skip until beggining of row is found
			while (fscanf(f, "%c", &c) && c != '[') { /* skip */ } 

			read_block(f, system.pS[block], row);

			//remove terminator
			fscanf(f, "]\n"); 
		}   
	}	
	
	return system;
}

/// serialize system
int write_mrhs_variable(FILE *f, MRHS_system system)
{
	int row, block, nrows, maxrhs;

	if (system.nblocks == 0) return 0;

	//number of lines for M
	nrows = system.pM[0].nrows;

	//print header
	fprintf(f, "%i %i\n", nrows, system.nblocks);
	for (block = 0; block < system.nblocks; block++)
	{
		fprintf(f, "%i %i\n", system.pS[block].nrows, system.pS[block].ncols);
	}		
	
	//print matrix
	for (row = 0; row < nrows; row++)
	{     
		fprintf(f, "[ ");
		for (block = 0; block < system.nblocks; block++)
		{
			print_block(f, system.pM[block], row);
			fprintf(f, " ");
		}
		fprintf(f, "]\n");
	}   

	// print RHS part
	for (block = 0; block < nblocks; block++)
	{
		fprintf(f, "\n");
		for (row = 0; row < maxrhs; row++)
		{     
			fprintf(f, "[");
			print_block(f, system.pS[block], row);
			fprintf(f, "]\n");
		}
	} 	
}

///print block matrix in user readable form
int print_mrhs(FILE *f, MRHS_system system)
{
	int row, block, nrows, maxrhs;

	if (system.nblocks == 0) return 0;

	//number of lines for M
	nrows = system.pM[0].nrows;

	//max number of lines for S
	maxrhs = system.pS[0].nrows;
	for (block = 1; block < system.nblocks; block++)
	{
		if (system.pS[block].nrows > maxrhs)
		{
			maxrhs = system.pS[block].nrows;
		}
	}
	
	//print matrix
	for (row = 0; row < nrows; row++)
	{     
		for (block = 0; block < system.nblocks; block++)
		{
			print_block(f, system.pM[block], row);
			fprintf(f, " ");
		}
		fprintf(f, "\n");
	}   

	//separator
	for (block = 0; block < system.nblocks; block++)
	{
		for (bit = 0; bit < system.pM[block].ncols; bit++)
		{
			fprintf(f, "-");
		}
		fprintf(f, " ");
	}
	fprintf(f, "\n");
	
	// print RHS part
	for (row = 0; row < maxrhs; row++)
	{     
		for (block = 0; block < nblocks; block++)
		{
			// no more rhs in given S
			if (i >= system.pS[block].nrows)
			{
				fprintf(f, "%*c", system.pS[block].ncols+1, ' ');
				continue;
			}   

			print_block(f, system.pS[block], row);
			fprintf(f, " ");
		}
		fprintf(f, "\n");
	}   
}

