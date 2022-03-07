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

#include "mrhs.hillc.h"

////////////////////////////////////////////////////////////////////////////////
// Macros and control structures

//0 -> stats only
//1 -> pretty print of stats
//2 -> include solutions
//3 -> include original system
//4 -> include statistics

#ifndef _VERBOSITY
 #define _VERBOSITY 2
#endif

#define MAXS 30

//setup of experiment and stats
typedef struct {
  int n,    // number of variables,      CMD LINE -n
      m,    // number of MRHS equations, CMD LINE -m
      l,    // dimension of RHS,         CMD LINE -l
      k;    // number of vectors in RHS, CMD LINE -k
  int seed; // seed for system generator, CMD LINE -s
  int seed2;// seed for solver, CMD LINE -S
  int rank; // row rank of the system (typically: rank == n)
  
  int d;    //system density
  
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


////////////////////////////////////////////////////////////////////////////////
// Utility functions

///fill in with random values
void random_bbm(_bbm *pbbm)
{
     int i, block;
     _block mask;
     
     for (block = 0; block < pbbm->nblocks; block++)
     {
       mask = (1u << pbbm->blocksizes[block]) - 1;
       for (i = 0; i < pbbm->nrows; i++)
       {     
           pbbm->rows[i][block] = rand() & mask;
       }
     }   
}

///fill in with random values, density = extra ones to place
void random_sparse_bbm(_bbm *pbbm, int density)
{
     int i, block, blocksize, row, col;
     _block mask;
     
     for (block = 0; block < pbbm->nblocks; block++)
     {
       //add random variables active in each column of a block
       blocksize = pbbm->blocksizes[block];
       mask = ONE;
       
       for (i = 0; i < blocksize; i++)
       {   
       	   row = rand() % pbbm->nrows;
       	   if (pbbm->rows[row][block] != ZERO)
       	   {
       	   	   //already set to another one
       	   	   i--; continue;
       	   }
       	   //set active variable
           pbbm->rows[row][block] = mask;
           mask <<= 1;
       } 
     }   
     for (i = 0; i < density; i++)
     {
       	   row = rand() % pbbm->nrows;
     	   block = rand() % pbbm->nblocks;
     	   col = rand() % pbbm->blocksizes[block];
           pbbm->rows[row][block] |= (ONE << col);     	   
     }
     
}

///print block matrix
void print(FILE* f, _bbm *pbbm, int pivots)
{
     int i, j, block;
     _block data;
     
     for (i = 0; i < pbbm->nrows; i++)
     {     
       for (block = 0; block < pbbm->nblocks; block++)
       {
          data = pbbm->rows[i][block]; 
          for (j = 0; j < pbbm->blocksizes[block]; j++, data>>=1)
          {
              fprintf(f, "%li", data&1);
          }
          fprintf(f, " ");
       }
       fprintf(f, "\n");
     }   
    
     if (!pivots) 
        return;

     for (i = 0; i < pbbm->nblocks; i++)
          printf("%*i%c", (int)pbbm->blocksizes[i], pbbm->pivots[i], 
                         i < pbbm->nblocks - 1 ? ' ' : '\n');
}

///print block matrix
void print_rhs(FILE* f, _bbm *prhs[], int nblocks)
{
     int i, j, block, maxrows;
     _block data;
     
     maxrows = prhs[0]->nrows;
     for (j = 1; j < nblocks; j++)
        if (prhs[j]->nrows > maxrows)
            maxrows = prhs[j]->nrows;
     
     for (i = 0; i < maxrows; i++)
     {     
       for (block = 0; block < nblocks; block++)
       {
          if (i >= prhs[block]->nrows)
          {
              fprintf(f, "%*c", (int)prhs[block]->blocksizes[0]+1, ' ');
              continue;
          }   
          
          data = prhs[block]->rows[i][0]; 
          for (j = 0; j < prhs[block]->blocksizes[0]; j++, data>>=1)
          {
              fprintf(f, "%li", data&1);
          }
          fprintf(f, " ");
       }
       fprintf(f, "\n");
     }   
}

///read system in sage format, first M by lines, then solutions by sets Si
int read_system_sage(FILE *f, _bbm *pbbm, _bbm *prhs[])
{
     int i, bit, block, read;
     
     //TODO: check correctness
     
     //list of rows of M
     for (i = 0; i < pbbm->nrows; i++)
     {   
         fscanf(f, "["); 
         for (block = 0; block < pbbm->nblocks; block++)
         {
             for (bit = 0; bit < pbbm->blocksizes[block]; bit++)
             {
                 fscanf(f, "%i", &read); 
                 pbbm->rows[i][block] ^= (((_block)read)<<bit);     //fixed block type
             }
         }   
         fscanf(f, "]\n"); 
     }
     
     //list of solutions S
     for (block = 0; block < pbbm->nblocks; block++)
     {
          for (i = 0; i < prhs[block]->nrows; i++)
          {     
             fscanf(f, "["); 
             for (bit = 0; bit < prhs[block]->blocksizes[0]; bit++)
             {
                 fscanf(f, "%i", &read); 
                 prhs[block]->rows[i][0] ^= (((_block)read)<<bit);        //fixed block type           
             }                
             fscanf(f, "]\n"); 
          }
     }   
     
     return 1;
}

///read system in sage format, first metadata, then M by lines, then solutions by sets Si
typedef struct {_bbm *pbbm; _bbm **prhs; } MRHS_system;

MRHS_system read_system_sage_new(FILE *f)
{
     MRHS_system system;
     char c;
     int i, bit, block, read, n, m, *k, suml, sumk;
     int *l;
     
     fscanf(f, "%i", &n);
     fscanf(f, "%i", &m);
    
     k = (int*)calloc(m, sizeof(int));
     l = (int*)calloc(m, sizeof(int));
     suml = 0; sumk = 0;
     for (i = 0; i < m; i++)
     {
         fscanf(f, "%i", &l[i]);
         fscanf(f, "%i", &k[i]);
         suml+=l[i];
         sumk+=k[i];
     }

     if (gp_experiment != NULL)
     {
        gp_experiment->n = n;
        gp_experiment->m = m;
        gp_experiment->l = suml/m;
        gp_experiment->k = sumk/m;
     }

     system.pbbm = create_bbm_new(n, m, l);
     system.prhs = (_bbm**) calloc(m, sizeof(_bbm*));
     //list of rows of M
     for (i = 0; i < system.pbbm->nrows; i++)
     {   
         while (fscanf(f, "%c", &c) && c != '[')
             ; 
         for (block = 0; block < system.pbbm->nblocks; block++)
         {
             for (bit = 0; bit < system.pbbm->blocksizes[block]; bit++)
             {
                 fscanf(f, "%1i", &read); 
                 system.pbbm->rows[i][block] ^= (((_block)read)<<bit);  //fixed block type
             }
         }   
         fscanf(f, "]\n"); 
     }
     
     //list of solutions S
     for (block = 0; block < system.pbbm->nblocks; block++)
     {
          system.prhs[block] = create_bbm(k[block], 1, l[block]);
          for (i = 0; i < system.prhs[block]->nrows; i++)
          {     
             while (fscanf(f, "%c", &c) && c != '[')
                ; 
             for (bit = 0; bit < system.prhs[block]->blocksizes[0]; bit++)
             {
                 fscanf(f, "%1i", &read); 
                 system.prhs[block]->rows[i][0] ^= (((_block)read)<<bit);     //fixed block type              
             }                
             fscanf(f, "]\n"); 
          }
     }   
     
     free(l);
     free(k);
     
     return system;
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
    fprintf(stderr, "\nUsage: %s [-n N] [-m M] [-l L] [-k K] [-s SEED] [-f FILE] [-e] [-t maxtime]\n", fn);
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
    setup->seed2 = -1; 
    setup->t = MAXS;
    setup->d = -1;
    
  while ((c = getopt (argc, argv, "ehk:l:m:n:s:S:f:t:d:")) != -1)
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
      case 'S':
        sscanf(optarg, "%i", &(setup->seed2));
        break;
      case 't':
        sscanf(optarg, "%lf", &(setup->t));
        break;
      case 'f':
        *pfn = optarg;
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
}



int main(int argc, char* argv[])
{
    //variables, equations, equation "degree", num rhs
    _experiment experiment;
    
    clock_t start, end;  
    _bbm *pbbm, **prhs;       
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
            //prhs[block] = create_bbm(experiment.k-1+rand()%3, 1, experiment.l);
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
    
    if (experiment.seed2 == -1)
        experiment.seed2 = time(0);
    srand(experiment.seed2);
    
#if (_VERBOSITY > 0)
    printf("Experimental setup, SEED = %08x, SEED2 = %08x \n", experiment.seed, experiment.seed2);
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
    
    if (!estimate)
    {   
        start = clock(); 
        //xors -> count of eval, total -> number of restarts
        experiment.count = solve(pbbm, prhs, start+experiment.t*CLOCKS_PER_SEC, &experiment.xors, &experiment.total);
        end = clock();
        experiment.t= (end-start)/(double)CLOCKS_PER_SEC;
    }

#if (_VERBOSITY > 0)
    printf("\nXORs: %lld Expected: %.0lf - %.0lf\n", 
               experiment.xors, experiment.xor2, experiment.xor1);
    printf("\nSolutions: %lld\nSearched %lld in %.3lf s, %e per sec\n", 
               experiment.count, experiment.total, experiment.t, experiment.total/experiment.t);
#endif
       
#if (_VERBOSITY == 0)

    //     SEED/SEED2   n   m   l  k rank count total time expected
    printf("%08x/%08x\t%i\t%i\t%i\t%i\t", 
		experiment.seed,experiment.seed2,experiment.n, experiment.m, experiment.l, experiment.k);
    printf("%i\t", 
		experiment.rank);
    printf("%lld\t%lld\t%lf\t%.0lf\t", 	
               experiment.count, experiment.total, experiment.t, experiment.expected);
    printf("%lld\t%.0lf\t%.0lf\n", 
               experiment.xors, experiment.xor2, experiment.xor1);
#endif
        
    for (block = 0; block < pbbm->nblocks; block++)
        free_bbm(prhs[block]); 
    free(prhs);
    free_bbm(pbbm);
     
    //system("pause");
    return 0;
}

//TODO: check solution, store all solutions?

