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
 #define _VERBOSITY 4
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
int parse_cmd(int argc, char *argv[], _experiment *setup, char **pfn, char **pfo);

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

int parse_cmd(int argc, char *argv[], _experiment *setup, char **pfn, char **pfo)
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
    
  while ((c = getopt (argc, argv, "ehk:l:m:n:s:S:f:o:t:d:")) != -1)
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
      case 'o':
        *pfo = optarg;
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
    char *fname = NULL;
    char *oname = NULL;
    FILE *f = NULL;
    FILE *fout = NULL;
    int estimate, block;
    MRHS_system system;

    estimate = parse_cmd(argc, argv, &experiment, &fname, &oname);
    
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
    if (oname != NULL)
    {
         fout = fopen(oname, "w");     
         if (fout == NULL)
         {
               fprintf(stderr, "Invalid file name: %s\n", oname);
               return -1;
         }
    }	
    //init stats
    init_stats(&experiment); 

    
    srand(experiment.seed);


    if (f == NULL)
    {
        system = create_mrhs_fixed(experiment.n, experiment.m, experiment.l, experiment.k);

		if (experiment.d == -1) 
		{
        	fill_mrhs_random(&system);
        }
        else
        {
        	fill_mrhs_random_sparse_extra(&system, experiment.d);
        }
    }
    else
    {
        system = read_mrhs_variable(f);
        fclose(f);
        
        experiment.m = system.nblocks;
        experiment.n = system.nblocks == 0 ? 0 : system.pM[0].nrows;
        experiment.l = system.nblocks == 0 ? 0 : system.pS[0].ncols;
        experiment.k = system.nblocks == 0 ? 0 : system.pS[0].nrows;
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
    print_mrhs(stdout, system);
#endif    

	if (fout != NULL)
	{
		write_mrhs_variable(fout, system);
		fclose(fout);
		fout = NULL;
	}
    
    if (!estimate)
    {   
        start = clock(); 
        //xors -> count of eval, total -> number of restarts
        experiment.count = solve(system, start+experiment.t*CLOCKS_PER_SEC, &experiment.xors, &experiment.total);
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
        
    clear_MRHS(&system);
     
    //system("pause");
    return 0;
}

//TODO: check solution, store all solutions?

