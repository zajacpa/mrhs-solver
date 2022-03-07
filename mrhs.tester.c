/**********************************
 * MRHS based solver
 * (C) 2016 Pavol Zajac  
 * 
 * I/O + testing module
 *
 * v1.3: refactoring + solution reporting + variable sized blocks
 * v1.4: minor fixes (reporting)
 * v1.6: fix input type - int to 64 bit when constructing matrix
 *  
 * Compilation: gcc -o mrhs mrhs.1.4.c tester.c -lm
 **********************************/
 
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "mrhs.solver.h"

////////////////////////////////////////////////////////////////////////////////
// Macros and control structures

//0 -> stats only
//1 -> pretty print of stats
//2 -> include solutions
//3 -> include original system
//4 -> include modified system
//5 -> include LUT
//6 -> include S*M
#ifndef _VERBOSITY
 #define _VERBOSITY 2
#endif

//setup of experiment and stats
typedef struct {
  int n,    // number of variables,      CMD LINE -n
      m,    // number of MRHS equations, CMD LINE -m
      l,    // dimension of RHS,         CMD LINE -l
      k;    // number of vectors in RHS, CMD LINE -k
  int seed; // seed for random generator, CMD LINE -s
  int rank; // row rank of the system (typically: rank == n)
  
  int d;    // system density 
  
  // number of solutions
  long long int count;
  // number of lookups (loops), predicted number of lookups
  long long int total; double expected;
  // counted number of xors, predicted number of xors (all), predicted number without adding zero-rows 
  long long int xors; double xor1, xor2;
  // measured time in seconds  
  double t;
} _experiment;

// Global pointer to experimental setup and stat reporting
static _experiment *gp_experiment = NULL;
// Sets gp_experiment to point to setup, and resets stats (count, total, xors, t)
void init_stats(_experiment *setup);

// Fills in experimental setup from command line arguments
//   Additionally: returns file name conatining eq. system, CMD LINE -f
//                 returns flag, whether run full experiment or just estimate, CMD LINE -e   
int parse_cmd(int argc, char *argv[], _experiment *setup, char **pfn);


///print block matrix solution - final solutions are on diagonal

void print_sol(FILE* f, _bbm *pbbm, ActiveListEntry* ale)
{
     int j, block, bitoffset;
     _block value;
     
     bitoffset = 0;
       for (block = 0; block < pbbm->nblocks; block++)
       {
          //nr = ale[block].u; 
          //value = (nr[bitoffset/MAXBLOCKSIZE]>>(bitoffset%MAXBLOCKSIZE))^ (nr[bitoffset/MAXBLOCKSIZE+1]<<(MAXBLOCKSIZE-bitoffset%MAXBLOCKSIZE));
          value = ale[block].val;
          
          for (j = 0; j < pbbm->blocksizes[block]; j++, value>>=1)
          {
              fprintf(f, "%lx", value&1);
          }
          fprintf(f, " ");
          bitoffset += pbbm->blocksizes[block];
       }
       fprintf(f, "\n");
    
}

void report_solution(long long int counter, _bbm *pbbm, ActiveListEntry* ale)
{
     printf("solution: %lli\n", counter);
     print_sol(stdout, pbbm, ale);     
}

_bbm *GlobalA = NULL;

//TODO: create function in solver to get solution y, and to multiply y*A
void report_solution_extract_y(long long int counter, _bbm *pbbm, ActiveListEntry* ale)
{
     int block, pivot, i;
     _block value, y[pbbm->nrows];
     
     printf("solution: %lli\n", counter);
     print_sol(stdout, pbbm, ale);   

     //TODO: create better functions for this...
     i = 0;
     printf("Vector y: ");     
     for (block = 0; block < pbbm->nblocks; block++)
     {
         value = ale[block].val;
         value>>=(pbbm->blocksizes[block]-pbbm->pivots[block]);
         for (pivot = 0; pivot < pbbm->pivots[block]; pivot++, value>>=1)
         {
             y[i++] = value&1;
             fprintf(stdout, "%lx", value&1);
         }
     }
     printf("\n");     
     
     printf("Vector x: ");     
     for (block = 0; block < GlobalA->nblocks; block++)
     {
         value = 0;
         for (i = 0; i < GlobalA->nrows; i++)
         {
            value ^= y[i] & GlobalA->rows[i][block];
         }
         fprintf(stdout, "%lx", value&1);
     }
     printf("\n");     
    
}



///print raw block data in hex
void print_raw_blocks(FILE* f, _block *blocks, int count, int bl)
{
     int i, j;
       for (i = 0; i < count; i++)
       {
       	  for (j = 0; j < bl; j++)
			fprintf(f, "%0lx", (blocks[i]>>j)&1);
       }
       fprintf(f, "\n");
}

///print block matrix solution - final solutions are on diagonal
void print_luts(FILE* f, ActiveListEntry* ale, int count, int blockcount, int bl)
{
     _block value;
     int block;
     _block size;
     TableEntry *p;
     
     for (block = 0; block < count; block++)
     {
        size = ale[block].mask + 1;
        fprintf(f, "LUT %i of size %lu:\n", block, size);
        
        for (value = 0; value < size; value++)
        {
        	p = ale[block].LUT[value];
            fprintf(f, "\t%p\n", p);
            while (p != NULL)
            {
                fprintf(f, "\t\t%lx\n", p->value);
                if (p->sm_row != NULL) print_raw_blocks(f, p->sm_row, blockcount, bl);
                p = p->next;
			}
            fprintf(f, "\n");               
        }
     }
     fprintf(f, "\n");   
}

void fill_random_rhs(_bbm *prhs[], int nblocks)
{
   int k, l, block, potential, r, count;
   _block x;
   
   for (block = 0; block < nblocks; block++)
   {
       l = prhs[block]->blocksizes[0];
       potential = (1ul)<<l;
       k = prhs[block]->nrows;
       count = 0;
       
       for (x = 0; potential > 0 && k > 0; x++, potential--)
       {
           //add to block?
           //  k -> need to add this amount, 
           r = rand() % potential;
           if (r < k)
           {
               prhs[block]->rows[count][0] = x;
               count++;
               k--; 
           }
       }
       
   }
}

void help(char* fn)
{
    fprintf(stderr, "\nUsage: %s [-n N] [-m M] [-l L] [-k K] [-s SEED] [-f FILE] [-e]\n", fn);
    fprintf(stderr, "   N = number of variables (def. 10)\n");    
    fprintf(stderr, "   M = number of MRHS eqs  (def. 10)\n");    
    fprintf(stderr, "   L = dimension of RHSs   (def. 3)\n");    
    fprintf(stderr, "   K = num. vectors in RHS (def. 4)\n\n");    

    fprintf(stderr, "FILE = file containing MRHS system \n      (if none, system is randomly generated using SEED)\n");    
    fprintf(stderr, "File format: METADATA {numbers N M L1 K1 .. Lm Km} \n");    
    fprintf(stderr, "           N  VECTORS of size M*SUM(Li) {rows of joint system matrix}\n");    
    fprintf(stderr, "           K1 VECTORS of size L1   {vectors in 1st RHS} \n");
    fprintf(stderr, "           K2 VECTORS of size L2   {vectors in 2nd RHS} \n");
    fprintf(stderr, "           ... \n");
    fprintf(stderr, "           Km VECTORS of size Lm   {vectors in m-th RHS} \n");
    fprintf(stderr, "        example VECTOR = [0 1 0 1 1 0] (size 6)\n\n");            

    fprintf(stderr, "[-e] = if enabled, SW only estimates complexity, does not solve the system\n\n");    
}

int parse_cmd(int argc, char *argv[], _experiment *setup, char **pfn)
{
   int c;
   opterr = 0;
   int estimate = 0;

    //default settings    
    setup->n = 10; 
    setup->m = 10; 
    setup->l = 3; 
    setup->k = 4; 
    setup->seed = -1; 
    setup->d = -1;

  while ((c = getopt (argc, argv, "ehk:l:m:n:s:f:t:d:")) != -1)
    switch (c)
      {
      case 'k':
        sscanf(optarg, "%i", &(setup->k));
        break;
      case 'l':
        sscanf(optarg, "%i", &(setup->l));
        break;
      case 'm':
        sscanf(optarg, "%i", &(setup->m));
        break;
      case 'n':
        sscanf(optarg, "%i", &(setup->n));
        break;
      case 'd':
        sscanf(optarg, "%i", &(setup->d));
        break;
      case 's':
        sscanf(optarg, "%i", &(setup->seed));
        break;
      case 'f':
        *pfn = optarg;
        break;
      case 't':
        //ignore, for compatibility
        break;
      case 'e':
        estimate = 1;
        break;
      case '?':
      case 'h':
        help(argv[0]);  
        exit(1);
      default:
        abort ();
      }


   return estimate;
}

//init stats
void init_stats(_experiment *setup)
{
    gp_experiment = setup;
    setup->xors  = 0; 
    setup->count = 0;
    setup->total = 0;
    setup->t     = 0.0;
}

int main(int argc, char* argv[])
{
    //variables, equations, equation "degree", num rhs
    _experiment experiment;
    ActiveListEntry* pActiveList;
    
    clock_t start, end;  
    _bbm *pbbm, **prhs, *pA = NULL;       
    char *fname = NULL;
    FILE *f = NULL;
    int estimate, block;
    MRHS_system system;

    estimate = parse_cmd(argc, argv, &experiment, &fname);
    
    if (experiment.seed == -1)
        experiment.seed = time(0);
    if (fname != NULL)
    {
         f = fopen(fname, "r");     
         if (f == NULL)
         {
               fprintf(stderr, "Invalid file name: %s\n", fname);
               return -1;
         }
    }
	
    //init stats
    init_stats(&experiment); 

    
    srand(experiment.seed);


    if (f == NULL)
    {
        pbbm = create_bbm(experiment.n, experiment.m, experiment.l);
        prhs = (_bbm**) calloc(pbbm->nblocks, sizeof(_bbm*));
        for (block = 0; block < pbbm->nblocks; block++)
            prhs[block] = create_bbm(experiment.k, 1, experiment.l);

		if (experiment.d == -1)
		{
        	random_bbm(pbbm);
        }
        else
        {
        	random_sparse_bbm(pbbm, experiment.d);
        }
        fill_random_rhs(prhs, pbbm->nblocks);
    }
    else
    {
        system = read_system_sage_new(f);
        pbbm = system.pbbm;
        prhs = system.prhs;
        fclose(f);
    }
    
#if (_VERBOSITY > 0)
    printf("Experimental setup, SEED = %08x \n", experiment.seed);
    printf("Num variables n = %i \n", experiment.n);
    printf("Num equations m = %i \n", experiment.m);
    printf("Block length  l = %i \n", experiment.l);
    printf("Block size    k = %i \n", experiment.k);
    if (f != NULL)
        printf("File: %s\n", fname);
#endif    

    
#if (_VERBOSITY > 2)
    print(stdout, pbbm, 0);
    fprintf(stdout,"-------------------------\n");
    print_rhs(stdout, prhs, pbbm->nblocks);
#endif    
    
#if (_VERBOSITY > 1)
    experiment.rank = echelonize(pbbm, prhs, &pA);
    GlobalA = pA;
 #if (_VERBOSITY > 3)
    fprintf(stdout,"Matrix A:\n");
    print(stdout, pA, 0);
 #endif
#else			 
    experiment.rank = echelonize(pbbm, prhs, NULL);
#endif
    pActiveList = prepare(pbbm, prhs);
    experiment.expected = get_expected(pbbm,prhs);

    experiment.xor1 = get_xor1(pbbm,prhs);
    experiment.xor2 = get_xor2(pbbm,prhs);

#if (_VERBOSITY > 0)
    printf("Expected count  = %.0lf \n", experiment.expected);
#endif    

#if (_VERBOSITY > 3)
    printf("\nRank: %i\n\n", experiment.rank);
    print(stdout, pbbm, 1);    
    print_rhs(stdout, prhs, pbbm->nblocks);
#endif

#if (_VERBOSITY > 4)
    printf("\n\nLUTS\n");
    print_luts(stdout, pActiveList, pbbm->nblocks, GET_NUM_BLOCKS_NEEDED(pbbm->ncols), MAXBLOCKSIZE);       
#endif
    
    if (!estimate)
    {   
        start = clock(); 
#if (_VERBOSITY > 1)
        experiment.total = solve(pActiveList, pbbm, &experiment.count, &experiment.xors, report_solution_extract_y);   
#else			 
        experiment.total = solve(pActiveList, pbbm, &experiment.count, &experiment.xors, NULL);   
#endif
        end = clock();
        experiment.t= (end-start)/(double)CLOCKS_PER_SEC;
    }
    else
    {
        experiment.total = experiment.expected;
        experiment.t = -1.0;
    }

    free_ales(pActiveList, pbbm->nblocks);

#if (_VERBOSITY > 0)
    printf("\nXORs: %lld Expected: %.0lf - %.0lf\n", 
               experiment.xors, experiment.xor2, experiment.xor1);
    printf("\nSolutions: %lld\nSearched %lld in %.3lf s, %e per sec\n", 
               experiment.count, experiment.total, experiment.t, experiment.total/experiment.t);
#endif
       
#if (_VERBOSITY == 0)

    //     SEED   n   m   l  k rank count total time expected
    printf("%08x\t%i\t%i\t%i\t%i\t", 
		experiment.seed,experiment.n, experiment.m, experiment.l, experiment.k);
    printf("%i\t", 
		experiment.rank);
    printf("%lld\t%lld\t%lf\t%.0lf\t", 	
               experiment.count, experiment.total, experiment.t, experiment.expected);
    printf("%lld\t%.0lf\t%.0lf\n", 
               experiment.xors, experiment.xor2, experiment.xor1);
#endif
        
    if (pA != NULL)
        free_bbm(pA);    
    for (block = 0; block < pbbm->nblocks; block++)
        free_bbm(prhs[block]); 
    free(prhs);
    free_bbm(pbbm);
     
        
    //system("pause");
    return 0;
}

//TODO: check solution, store all solutions?

