/**********************************
 * MRHS based solver
 * (C) 2016 Pavol Zajac  
 * 
 * library file
 *
 * v1.3: refactoring + solution reporting + variable sized blocks
 * v1.4: minor update: xor counting + valgrind / memfixes
 * v1.5: bugfix in swapcols
 * v1.6: swapcols again...
 **********************************/
 //TODO: create simple function: 
 //                input MRHS system, stat pointer, callback
 //                output: stats, 
 //                callback per each solution: input - solution, output: stop/continue

#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "mrhs.solver.h"

////////////////////////////////////////////////////////////////////////////////
// BBM Echelon form

///row to ^= to_add 
void add_row(_bbm *pbbm, int to, int to_add)
{
    int i; 
    for (i = 0; i < pbbm->nblocks; i++) {
        pbbm->rows[to][i] ^= pbbm->rows[to_add][i];
    }
}

///swaps rows i and j
void swap_row(_bbm *pbbm, int i, int j)
{
     _block* tmp   = pbbm->rows[i]; 
     pbbm->rows[i] = pbbm->rows[j];
     pbbm->rows[j] = tmp;
}

///swaps cols in a block
void swap_cols(_bbm *pbbm, int block, int col1, int col2)
{
     int row; 
     _block data;
     _block mask = ((ONE)<<col1) | ((ONE)<<col2);  //BUGFIX: changed from 1u to (_block)1 //constant ONE
     
     for (row = 0; row < pbbm->nrows; row++)
     {
         data  = pbbm->rows[row][block];
         //if both are same, no change, otherwise change both of them
         data ^= ( ((data >> col1)^(data>>col2))&ONE )*mask; //constant ONE
         
         pbbm->rows[row][block] = data;
     }
}

///find first non-zero element in given block/mask, from given row index 
/// -1 if no such exists
int find_pivot(_bbm *pbbm, int block, _block mask, int from)
{
    while (from < pbbm->nrows)
    {
        if (pbbm->rows[from][block] & mask)
            return from;
        from++;
    }
    return -1;
}

/// Creates an echelon form of matrix and computes number of pivots
///  swaps columns so that pivots form identity matrix at the start of block
/// if rhs pointer is non-zero, swaps columns also in rhs(must have same blocks)
/// NOTE: pivots are moved to MSB part, so that LSB part can be used as an index
//TODO: rhs is a list of RHSs, each with a different number of rows 
int echelonize(_bbm *pbbm, _bbm *rhs[], _bbm **pA)
{
     //ASSERT (rhs != NULL)
     //ASSERT: rhs is an array of matrices, single block each, variable number of rows
     //ASSERT: len(rhs) = pbbm->nblocks
     int rank = 0, block = 0, i, j, k, pivot, oldrank = 0, offset;
     _block mask, mask2;
     
     if (pA != NULL)
     {
        //create unit matrix "A" - for getting back the solution
        *pA = create_bbm(pbbm->nrows, pbbm->nrows, 1);
        for (block = 0; block < pbbm->nrows; block++)
            (*pA)->rows[block][block] = ONE;    
        
     }
     block = 0;
     while (rank < pbbm->nrows && block < pbbm->nblocks)
     {     
       //TODO: ASSERT(pbbm->pivots[block] == 0)
       pbbm->pivots[block] = 0;

       for (j = 0, mask = ONE; j < pbbm->blocksizes[block]; j++, mask <<= 1)
       {
           pivot = find_pivot(pbbm, block, mask, rank);
           
           if (pivot < 0)
           {
              continue;
           }
           
           //move pivot row to required position
           if (pivot != rank)
           {
               swap_row(pbbm, pivot, rank);
               if (pA != NULL)
                   swap_row(*pA, pivot, rank);
           }
               
           //now reduce with row at rank index
           for (i = 0; i < pbbm->nrows; i++)
           {
               if (i == rank) continue;
               if (pbbm->rows[i][block]&mask)
               {
                   add_row(pbbm, i, rank);
                   if (pA != NULL)
                       add_row(*pA, i, rank);
               }
           }
                     
           //perform swap if required
           if (j > pbbm->pivots[block])
           {
               swap_cols(pbbm, block, j, pbbm->pivots[block]);
               swap_cols(rhs[block], 0, j, pbbm->pivots[block]);
           }
           //update npivots
           pbbm->pivots[block]++;

           rank++;    
       }      

       //all done, zero out redundant part of block's matrix
       // we get (I | 0) in pivot part of the block
       for (pivot = 0, mask = ONE; //first pivot in LSB
                     pivot < pbbm->pivots[block]; pivot++, mask <<= 1)
       {
           for (j = pbbm->pivots[block]; j < pbbm->blocksizes[block]; j++)
           {
               mask2 = ONE << j;
               if ( (mask2&pbbm->rows[oldrank+pivot][block]) )
               {
                    //remove selected column from each sol
                    pbbm->rows[oldrank+pivot][block] ^= mask2;
                    for (k = 0; k < rhs[block]->nrows; k++)
                         rhs[block]->rows[k][0] ^= (mask&rhs[block]->rows[k][0])<<(j-pivot);
               }
           }
       }
   
       //move pivots to MSB part
       offset = pbbm->blocksizes[block] - pbbm->pivots[block];
       if (offset > 0)
       {
           for (j = pbbm->pivots[block] - 1; j >= 0; j--)
           {
               swap_cols(pbbm, block, j, j + offset);
               swap_cols(rhs[block], 0, j, j + offset);
           }
       }

                
       oldrank = rank;
       
       block++;
     }
    
     return rank;
}


/**************************************************************************
 * PRECOMP:
 *  1. Compute echelon form of M by row operations
 *  2. In each block: use column operations to bring last rows to form G = I|0, 
 *  3. For each s_i \in S, compute u_i = [s_i] M, where [s_i] is pivot part
 *
 ***************************************************************************/

//multiply part of bbm matrix - from offset, by coeff, adds result to out
// compress result into output block
// returns first non-zero block
int multiply_add(_block out[], _block coeff, _bbm *pbbm, int offset)
{
     int i, bitoffset, blocksize, ix, off;
     _block to_add, to_add_part2;
     
     for ( ; coeff && offset < pbbm->nrows; coeff >>= 1, offset++)
     {
         if ((coeff & ONE) == 0)
             continue;
         
	 bitoffset = 0;    
         for (i = 0; i < pbbm->nblocks; i++)
         {
             to_add = pbbm->rows[offset][i];
             
             blocksize = pbbm->blocksizes[i];
             ix = bitoffset / MAXBLOCKSIZE;   //array index
             off = bitoffset % MAXBLOCKSIZE;  //bit index in array
             
             out[ix] ^= ((to_add) << off);
             
             //if: to prevent write outside array
             to_add_part2 = (((to_add) >> (MAXBLOCKSIZE-1-off))>>1);
             if (to_add_part2 != 0)
                out[ix+1] ^= to_add_part2;	//could be out of array, but then we would be adding zero...
             
             bitoffset += blocksize;
         }
     }
     //use last remembered ix
     for (i = 0; i <= ix; i++)
	if (out[i] != 0)
           break;
     return i;      
}


//This function is used to disallow duplicate entries in LUTs
//TODO: add compile time define to turn this off
int contains(_block value, TableEntry *LUT){
	TableEntry *current = LUT;
	while(current){
		if(value == current->value)
			return 1;
		current = current->next;
	}
	return 0;
}

//PRE: pbbm and prhs prepared by echelonize
//TODO?: variable block sizes - this should already work 
//WORKAROUND: allows variable number of rhs by removing duplicate entries
ActiveListEntry* prepare(_bbm *pbbm, _bbm *prhs[])
{
    int block, rhs, offset, k, r;
    _block size, index, value;
    ActiveListEntry *pList;
    int blocklen = GET_NUM_BLOCKS_NEEDED(pbbm->ncols);  //pbbm->nblocks; //(*pbbm->blocksizes[0]/MAXBLOCKSIZE);
    
    //PRE: pbbm has echelon form, with (I|0) in blocks with free pivots
    //      pivots are stored from LSB bits
    //      prhs are sorted, pivot part is in upper (MSB) part - values
    //                       upper part corresponds to (0 target) parities  - LUT keys
    
    //allocate list of entries for search algorithm
    pList = (ActiveListEntry*) calloc(pbbm->nblocks, sizeof(ActiveListEntry));
      
    // this part precomputes list's lookup tables 
    //   along with S * M 
    offset = 0; k = 0;
    for (block = 0; block < pbbm->nblocks; block++)
    {
    //TODO: beware, dangerous alloc! can reach 2^blocksize
        //index mask ... 
        r = pbbm->blocksizes[block] - pbbm->pivots[block];
        size = ONE << r;     //size of LUT
        pList[block].mask = (size - 1);  //if r == 0 -> 0, else r ones
        pList[block].LUT = (TableEntry**) calloc(size, sizeof(TableEntry*));
        //pList[block].next = NULL;    //calloc...
        //pList[block].sol_ix = 0;
                
        for (rhs = 0; rhs < prhs[block]->nrows; rhs++)
        {
            value = prhs[block]->rows[rhs][0];
            index = value & (pList[block].mask);
            value ^= index;                       //remove lower part
            
            //check whether the value is in the lut list
            //prevents duplicates
            if(contains(value, pList[block].LUT[index]) == 1)
		continue;

            //TODO: allow more flexibility, including some sort order in LUT 
            
            //create new entry in linked list
            TableEntry* nte = (TableEntry*) calloc(1, sizeof(TableEntry));
            nte->value = value;
            nte->next = pList[block].LUT[index];
            pList[block].LUT[index] = nte;           
            
            if (value != 0)
            {
                //zero-out corresponding rows of S*M
                //memset(psm->rows[k], 0, sizeof(_block)*prhs->nblocks);
            	nte->sm_row = (_block*) calloc(blocklen, sizeof(_block));
                //if there are free pivots, compute corresponding s_i * M
                nte->first = multiply_add(nte->sm_row, 
                                  (value)>>r, //move it back
                                  pbbm, offset);
                                  
                k++;
            }
            else
            {
            	nte->sm_row = NULL;
			}
        }
        offset += pbbm->pivots[block];
    }
    
    return pList;
}

void free_ales(ActiveListEntry* ale, int count)
{
     int i;
     _block j;
     TableEntry* next, *nextnext;
     for (i = 0; i < count; i++)
     {
         for (j = 0; j <= ale[i].mask; j++)
         {
             next = ale[i].LUT[j];
             while (next != NULL)
             {
                   nextnext = next->next;
                   if (next->sm_row != NULL)
                		free(next->sm_row);
                   free(next);
                   next = nextnext;
             }
         }
         free(ale[i].LUT);
     }
     free(ale);
}




//TODO: variable block sizes, variable number of rhs
long long int solve_it(ActiveListEntry* ale, _bbm *pbbm, int block, 
                       _block* sol_stack, long long int *pCount, long long int *pXors,
                        sol_rep_fn_t report_solution)
{
    long long int count = 0;
    long long int xors = 0;
    long long int total = 0;
    int b;
    _block index;
    TableEntry * active;
    _block *nr, *or, *ar, value;  //new row, old row, active row, value from u
    int blocklen = GET_NUM_BLOCKS_NEEDED(pbbm->ncols);  //pbbm->nblocks; //(*pbbm->blocksizes[0]/MAXBLOCKSIZE);
    //int blocksize = pbbm->blocksizes[0];  //TODO: variable block sizes...
    int bitoffset;
    
    nr = or = sol_stack;
    
    //ASSERT: block == 1, sol->pivots == 0 for each block
    bitoffset = 0;
    while (block >= 0)
    {
        active = ale[block].next;
        
        //no more to process
        if (active == NULL)
        {
             //backtrack
             block--; 
             
             //TODO: possible optimalization: allocate one block before blocksizes to allow negative index
             if (block >= 0)
                 bitoffset -= pbbm->blocksizes[block];
             continue;
        }

        //reporting
        ++total;
        
        //prepare stack for next solution
        ale[block].next = active->next;
        
        //are we at the end?
        if (block == pbbm->nblocks-1)
        {
             //solution found    
             ++count;
		
		 	//recompute solution
            ale[block].val = (ale[block].val& ale[block].mask)^active->value;
			 
            if (report_solution!=NULL)
                report_solution(count, pbbm, ale);

            //break;
	    //TODO: why does it generate infinite loop...
            continue;
        }

        //no change required in solution if active.value == 0, else add sm
        if (active->value == 0)
        {
            //recompute solution
            ale[block].val = (ale[block].val& ale[block].mask)^active->value;
            
            ale[block+1].u = ale[block].u;
            nr = ale[block].u;
        }
        else
        {
            or = ale[block].u;
            nr = or+blocklen;
            ale[block+1].u = nr;
            
            //add to previous solution, before "block" all zeroes
            ar = active->sm_row;
            
            //or[block] = (or[block]& ale[block].mask)^active->value;
            ale[block].val = (ale[block].val& ale[block].mask)^active->value;
            
            //add block to u
            for (b = active->first; b < blocklen; b++)
            {
                 nr[b] = or[b]^ar[b];
                 //reporting:
                 xors++;
            }
        }
        
        //get LUT index
        bitoffset += pbbm->blocksizes[block];
        block++;
        
        //get bits from u...
        //bitoffset = (block*blocksize);
        value = (nr[bitoffset/MAXBLOCKSIZE]>>(bitoffset%MAXBLOCKSIZE))^ ((nr[bitoffset/MAXBLOCKSIZE+1]<<(MAXBLOCKSIZE-1-bitoffset%MAXBLOCKSIZE))<<1);
        index = value & ale[block].mask;
        ale[block].next = ale[block].LUT[index];
    	ale[block].val = index;

        //gp_experiment->lookups++;
    }
    //all done, back to root
    if (pCount != NULL)
        *pCount += count;
    if (pXors != NULL)
        *pXors += xors;
    
    return total;
}

//front end to non-recursive call
//TODO: for multiprocessing, fork can be used and new process created for each rhs
//TODO: for threading, sol must be created for each rhs/thread 
long long int solve(ActiveListEntry* ale, _bbm *pbbm, long long int *pCount, long long int *pXors, sol_rep_fn_t report_solution)
{
    long long int total = 0;
    
    //TODO: variable block size
    int blocklen = GET_NUM_BLOCKS_NEEDED(pbbm->ncols); //pbbm->nblocks; //(*pbbm->blocksizes[0]/MAXBLOCKSIZE);
	_block* solstack = (_block*) calloc(pbbm->nblocks*blocklen, sizeof(_block));
   
    if (pCount != NULL)
        *pCount = 0;
    
    ale[0].u = solstack;
    ale[0].next = ale[0].LUT[0];
    ale[0].val = 0;
    
    //redundant - stored in total
    //gp_experiment->lookups++;

    total = solve_it(ale, pbbm, 0, solstack, pCount, pXors, report_solution);     
    //free(myword);
    free(solstack);
    
    return total;
}

#if (_VERBOSITY > 4)
#include <stdio.h>
#endif

///formula from article Ntotal
/// sum ( prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_expected(_bbm *pbbm, _bbm *prhs[])
{
    double sum = 0, prod;
    int i,j;
    for (i = 1; i < pbbm->nblocks; i++) 
    {
    	prod = 1;
    	for (j = 0; j < i; j++)
    	{
    		prod *= prhs[j]->nrows;
    		prod *= (ONE<<pbbm->pivots[j]);
    		prod /= (ONE<<pbbm->blocksizes[j]);    		
    	}
    	sum += prod;
#if (_VERBOSITY > 4)
    	fprintf(stderr, "EXP %i: %lf %lf\n", i, prod, sum);
#endif
    }
    return sum;	
}

///formula from article Nxor
/// sum ( (m-i+1) prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_xor1(_bbm *pbbm, _bbm *prhs[])
{
    double sum = 0, prod;
    int i,j;
    for (i = 1; i < pbbm->nblocks; i++) 
    {
    	prod = 1;
    	for (j = 0; j < i; j++)
    	{
    		prod *= prhs[j]->nrows;
    		prod *= (ONE<<pbbm->pivots[j]);
    		prod /= (ONE<<pbbm->blocksizes[j]);    		
    	}
    	sum += (ceil((pbbm->nblocks-i)*pbbm->blocksizes[i]/(double)MAXBLOCKSIZE) * prod);  
		//m-i+1 = m-(i-1), i is zero-indexed in C, recomputed to paralel block size processing
#if (_VERBOSITY > 4)
    	fprintf(stderr, "XOR1 %i: %lf %lf\n", i, prod, sum);
#endif
    }
    return sum;	
}

///formula from article Nxored
/// sum ( (1-2^(-p{i-1})) (m-i+1) prod(|S_j|*2^(pj-lj) j=1 to i-1)  i = 2 to m)
double get_xor2(_bbm *pbbm, _bbm *prhs[])
{
    double sum = 0, prod, c;
    int i,j;
    for (i = 1; i < pbbm->nblocks; i++) 
    {
    	prod = 1;
    	for (j = 0; j < i; j++)
    	{
    		prod *= prhs[j]->nrows;
    		prod *= (ONE<<pbbm->pivots[j]);
    		prod /= (ONE<<pbbm->blocksizes[j]);    		
    	}
        c = (ONE<<pbbm->pivots[i-1]);
        c = (c-1)/c;
    	sum += (ceil(c*(pbbm->nblocks-i)*pbbm->blocksizes[i]/(double)MAXBLOCKSIZE) * prod);  
		//m-i+1 = m-(i-1), i is zero-indexed in C, recomputed to paralel block size processing
#if (_VERBOSITY > 4)
    	fprintf(stderr, "XOR2 %i: %lf %lf\n", i, prod, sum);
#endif
    }
    return sum;	
}
