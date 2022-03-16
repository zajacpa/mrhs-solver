/**********************************
 * MRHS based solver
 * (C) 2016 Pavol Zajac
 *
 * I/O + testing module
 *
 * v1.3: refactoring + solution reporting + variable sized blocks
 * v1.4: minor fixes (reporting)
 * v1.6: fix input type - int to 64 bit when constructing matrix
 * v1.7: refactoring + integration of HC
 *
 * Compilation: make...
 **********************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "mrhs.bv.h"
#include "mrhs.hillc.h"
#include "mrhs.rz.h"


/// ////////////////////////////////////////////////////////////////////
/// Stats and reporting

//0 -> stats only
//1 -> pretty print of stats
//2 -> include solutions
//3 -> include original system
//4 -> include statistics

#ifndef _VERBOSITY
 #define _VERBOSITY 4
#endif

#define RZ_SOLVER_TYPE 1
#define HC_SOLVER_TYPE 2


typedef struct {
  int rank; // row rank of the system (typically: rank == n)

  // number of solutions
  long long int count;

  // number of lookups (loops),
  long long int total;

  // predicted number of lookups
  double expected;

  // counted number of xors
  long long int xors;

  // predicted number of xors (all), predicted number without adding zero-rows
  double xor1, xor2;

  // measured time in seconds
  double t;
} _stats;

// Global pointer to experimental setup and stat reporting
static _stats *gp_stats = NULL;
// Sets gp_experiment to point to setup, and resets stats (count, total, xors, t)
void init_stats(_stats *stats);

//init stats
void init_stats(_stats *stats)
{
	stats->rank     = 0;

    stats->xors     = 0;
    stats->count    = 0;
    stats->total    = 0;

	stats->expected = 0.0;
    stats->xor1     = 0.0;
    stats->xor2     = 0.0;
    stats->t        = 0.0;

    gp_stats = stats;
}

/// ////////////////////////////////////////////////////////////////////
/// Command line interface


//setup of experiments
typedef struct {
  int n,    // number of variables,      CMD LINE -n
      m,    // number of MRHS equations, CMD LINE -m
      l,    // dimension of RHS,         CMD LINE -l
      k;    // number of vectors in RHS, CMD LINE -k
  int seed; // seed for system generator, CMD LINE -s
  int seed2;// seed for solver, CMD LINE -S

  int d;       //system density
  double maxt; //time limit in seconds
  int solver;   //solver type
  int compress;  //enable compression of system

  char *in;    // system  input file
  char *out;   // system output file
  FILE *fsols; // open output file for solutions
} _experiment;

// Fills in experimental setup from command line arguments
//   Additionally: returns file name conatining eq. system, CMD LINE -f
//                 returns flag, whether run full experiment or just estimate, CMD LINE -e
int parse_cmd(int argc, char *argv[], _experiment *setup);

void help(char* fn)
{
    fprintf(stderr, "\nUsage: %s [-n N] [-m M] [-l L] [-k K] [-s SEED] [-S SED2] [-f FILE] [-o OUT] [-e TYPE] [-t MAXT]\n", fn);
    fprintf(stderr, "   N = number of variables (def. 10)\n");
    fprintf(stderr, "   M = number of MRHS eqs  (def. 10)\n");
    fprintf(stderr, "   L = dimension of RHSs   (def. 3)\n");
    fprintf(stderr, "   K = num. vectors in RHS (def. 4)\n\n");
    fprintf(stderr, "MAXT = time limit (in seconds)\n");
    fprintf(stderr, "SEED = randomness seed for MRHS system\n");
    fprintf(stderr, "SED2 = randomness seed for computation\n\n");

    fprintf(stderr, "TYPE = solver type: 0=no solver, %d=Raddum-Zajac, %d=HC\n\n", RZ_SOLVER_TYPE, HC_SOLVER_TYPE);

    fprintf(stderr, "FILE = file containing MRHS system \n      (if none, system is randomly generated using SEED)\n");
    fprintf(stderr, "OUT  = file to write out generated MRHS system \n\n");
    fprintf(stderr, "File format: METADATA {numbers N M L1 K1 .. Lm Km} \n");
    fprintf(stderr, "           N  VECTORS of size M*SUM(Li) {rows of joint system matrix}\n");
    fprintf(stderr, "           K1 VECTORS of size L1   {vectors in 1st RHS} \n");
    fprintf(stderr, "           K2 VECTORS of size L2   {vectors in 2nd RHS} \n");
    fprintf(stderr, "           ... \n");
    fprintf(stderr, "           Km VECTORS of size Lm   {vectors in m-th RHS} \n");
    fprintf(stderr, "        example VECTOR = [0 1 0 1 1 0] (size 6)\n\n");

    //fprintf(stderr, "[-e] = if enabled, SW only estimates complexity, does not solve the system\n\n");
}

void set_default_experiment(_experiment *setup)
{
    //default settings
    setup->n     = 10;   //system size
    setup->m     = 10;
    setup->l     = 3;
    setup->k     = 4;

    setup->seed  = -1;   //time based seeds
    setup->seed2 = -1;

    setup->maxt  = 1.0;  //1s max for HC
    setup->d     = -1;   //dense matrix
    setup->solver = RZ_SOLVER_TYPE;    //try to solve with RZ solver
    setup->compress = 0;  //no equation compression

    setup->in    = NULL; //no input/output
    setup->out   = NULL;
    setup->fsols = NULL; //TODO...
}

int parse_cmd(int argc, char *argv[], _experiment *setup)
{
   int c;
   opterr = 0;

   set_default_experiment(setup);

   while ((c = getopt (argc, argv, "ce:hk:l:m:n:s:S:f:o:t:d:")) != -1)
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
        sscanf(optarg, "%lf", &(setup->maxt));
        break;
      case 'f':
        setup->in = optarg;
        break;
      case 'o':
        setup->out = optarg;
        break;
      case 'e':
        sscanf(optarg, "%i", &(setup->solver));
        break;
      case 'c':
        setup->compress = 1;
        break;
      case '?':
      case 'h':
        help(argv[0]);
        exit(1);
      default:
        abort ();
      }

   return 1;
}

int prepare_system(MRHS_system *system, _experiment *setup)
{
    FILE *f = NULL;
    FILE *fout = NULL;

    //check I/O files
    if (setup->in != NULL)
    {
		//read MRHS system from file
        f = fopen(setup->in, "r");
        if (f == NULL)
        {
           fprintf(stderr, "Invalid file name: %s\n", setup->in);
           return 0;
        }

        *system = read_mrhs_variable(f);
        fclose(f);

        setup->m = system->nblocks;
        setup->n = system->nblocks == 0 ? 0 : system->pM[0].nrows;
        setup->l = system->nblocks == 0 ? 0 : system->pS[0].ncols;
        setup->k = system->nblocks == 0 ? 0 : system->pS[0].nrows;
    }
    else
    {
		//create random system
		if (setup->seed == -1)
			setup->seed = time(0);
		srand(setup->seed);

		//empty system:
        *system = create_mrhs_fixed(setup->n, setup->m, setup->l, setup->k);

		//dense or sparse?
		if (setup->d == -1)
		{
        	fill_mrhs_random(system);
        }
        else
        {
        	fill_mrhs_random_sparse_extra(system, setup->d);
        }
	}

    //report system ?
    if (setup->out != NULL)
    {
         fout = fopen(setup->out, "w");
         if (fout == NULL)
         {
               fprintf(stderr, "Invalid file name: %s\n", setup->out);
               return 0;
         }
         write_mrhs_variable(fout, *system);
         setup->fsols = fout;
    }

	return 1;
}


int main(int argc, char* argv[])
{
	//working with this system
    MRHS_system system;
    _bv *results = NULL;

    //input settings: variables, equations, equation "degree", num rhs
    _experiment experiment;

    //output statistics:
    _stats stats;

    //time and IO
    clock_t start, end;

	//prepare parameters
    if (!parse_cmd(argc, argv, &experiment))
		return -1;

 	if (!prepare_system(&system, &experiment))
		return -2;

    //init stats
    init_stats(&stats);

    //init random generator for experiments
    if (experiment.seed2 == -1)
        experiment.seed2 = time(0);
    srand(experiment.seed2);

#if (_VERBOSITY > 0)
    printf("Experimental setup, SEED = %08x, SEED2 = %08x \n", experiment.seed, experiment.seed2);
    if (experiment.in != NULL)
        printf("Input file: %s\n", experiment.in);
    printf("Num variables n = %i \n", experiment.n);
    printf("Num equations m = %i \n", experiment.m);
    printf("Block length  l = %i \n", experiment.l);
    printf("Block size    k = %i \n", experiment.k);
    if (experiment.out != NULL)
        printf("System stored to: %s\n", experiment.out);
#endif

#if (_VERBOSITY > 2)
    printf("\nInitial MRHS system: \n");
    print_mrhs(stdout, system);
    printf("\n");
#endif

    if (experiment.compress)
    {
        //make linear equation substitutions
        remove_linear(&system);
        remove_empty(&system);
    #if (_VERBOSITY > 2)
        printf("\nCompressed MRHS system: \n");
        print_mrhs(stdout, system);
        printf("\n");
    #endif
    }


	// run the experiment

	start = clock();
	//xors -> count of eval, total -> number of restarts
    switch (experiment.solver)
    {
	case HC_SOLVER_TYPE:
        stats.count = solve_hc(&system, &results, start+experiment.maxt*CLOCKS_PER_SEC, &stats.xors, &stats.total);
        break;
	case RZ_SOLVER_TYPE:
        stats.count = solve_rz(&system, &results, start+experiment.maxt*CLOCKS_PER_SEC, &stats.xors, &stats.total);
        break;
    }
	end = clock();
	stats.t= (end-start)/(double)CLOCKS_PER_SEC;

	// post processing: report results and clear data structures

	if (experiment.fsols != NULL && results != NULL)
	{
		for (int i = 0; i < stats.count; i++)
		{
			fprintf(experiment.fsols, "\nx ");
			print_bv(&results[i], experiment.fsols);
		}

	}
	if (experiment.fsols != NULL)
	{
		fclose(experiment.fsols);
	}

	if (results != NULL)
	{

		for (int i = 0; i < stats.count; i++)
		{
#if (_VERBOSITY > 1)
			fprintf(stdout, "\nSolution %i: ", i+1);
			print_bv(&results[i], stdout);
#endif
			clear_bv(&results[i]);
		}
#if (_VERBOSITY > 1)
		fprintf(stdout, "\n");
#endif

		free(results);
	}

	clear_MRHS(&system);

	// post processing, report statistics

#if (_VERBOSITY > 0)
    printf("\nXORs: %lld Expected: %.0lf - %.0lf\n",
               stats.xors, stats.xor2, stats.xor1);
    printf("\nSolutions: %lld\nSearched %lld in %.3lf s, %e per sec\n",
               stats.count, stats.total, stats.t, stats.total/stats.t);
#endif

#if (_VERBOSITY == 0)

    //     SEED/SEED2   n   m   l  k rank count total time expected
    printf("%08x/%08x\t%i\t%i\t%i\t%i\t",
		experiment.seed,experiment.seed2,experiment.n, experiment.m, experiment.l, experiment.k);
    printf("%i\t",
		stats.rank);
    printf("%lld\t%lld\t%lf\t%.0lf\t",
               stats.count, stats.total, stats.t, stats.expected);
    printf("%lld\t%.0lf\t%.0lf\n",
               stats.xors, stats.xor2, stats.xor1);
#endif

    //system("pause");
    return 0;
}

//TODO: check solution, store all solutions?
