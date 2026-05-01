
/* random number generator */
#define _MT


/* threshold for converting to/from integer */
#define EPS 1e-10

/*
  Include all dependencies
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Library that measures dates and times
#include <assert.h>
#include "mt19937-64.c" // random number generator
#include "util.c" // utilities
#include "spatial_hashing.c" // spatial hashing
#include "inhom_fluct.c" // inhomogeneous fluctuation landscape


// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param {
  double Lx;                  // system size in x
  double Ly;                  // system size in y
  long N;                     // integers can be int, long, long int.
  double dt;                  // Time step
  double tf;                  // end of simulation (real time, not # steps)
  long long seed;             // Seed of random number generator

  fparam_cub TParam;          // parameters of temperature landscape (included in inhom_fluct.c)
  double sqrt2dt;
  char* Ts;
  char* Txs;

  double a;                   // interaction length
  double rmax;                // interaction cutoff
  double k;                   // interaction strength
  double rmax2;               // interaction cut-off (squared)
  long NxBox;                 // number of boxes along x
  long NyBox;                 // number of boxes along y
  double c1,c2,c3;            // force parameters
  box ***NeighbouringBoxes;

  double UpdateInterProg;     // Interval between 2 updates of progress bar
  double StoreInterPos;       // Interval between 2 storages of the position
  double binwidthx;           // width of position profile bins (only on x axis)
  double binwidthy;           // width of position profile bins (only on y axis)
  int Nbinx;                  // number of bins for position profiles (=Lx/binwidth)
  int Nbiny;                  // number of bins for position profiles (=Ly/binwidth)
  double StoreInterProf;      // Interval between 2 storages of the pos
                              // (and, momentum, if underdamped) profile
  int NstepProf;              // number of steps to use when storing histogram
  double UpdateInterProf;     // Interval between 2 updates of profile

  FILE *output_param;         // File where parameters are stored
  FILE *output_pos;           // File where position is stored
  FILE *output_Tprof;         // File where temperature profile is stored
  FILE* output_sigmaIKprof;   // file where profile of stress tensor is stored
  FILE* output_EPR;
} param;


// Structure particles which contains all the data needed to characterize the
// state of a particle
typedef struct particle
{
  double x;
  double y;     // position of particles
  double dx;    // x increment
  double dy;    // y increment
  double T;     // temperature
  double ifx;   // interaction force in x direction
  double ify;   // interaction force in y direction
  int bi;       // box membership
  int bj;       // box membership
  double* ifxs; // list of x interaction forces due to all others
  double* ifys; // list of y interaction forces due to all others
} particle;


typedef struct state{
  clock_t start_time;              // real time of simulation start
  double _time;                    // Current simulation time
  double prev_percentage;          // previous simulation progress
  double NextUpdateProg;           // next time to print simulation progress

  double NextStorePos;             // Next time after which the position is stored
  double NextStoreProf;            // Next time after which the position profile is stored
  double NextUpdateProf;           // Next time to add to density/magnetization profile
} state;


typedef struct data{
  double*** sigmaIKprof;           // sigmaIKprof[i][j][k] = stress component k at xi,xj
  double** Tprof;                  // Tprof[i][j] = rho*T at xi,xj
  double** EPR;                    // EPR[i][j] = EPR at xi,xj
} data;

#include "PFs.c"
#include "./OD-PBPs-2D-functions.c" // our methods

int main (int argc, char *argv[]) {

  /* declaration of variables */
  particle *Particles; // Array that contains the state of all particles
  param Param;         // Param structure that contains information given to functions
  state State;         // state of simulation
  data Data;           // data arrays

  long **Boxes;
  long *Neighbours;

  /* initializing all variables */
  Initialize_parameters(argc, argv, &Param, &State, &Data, &Particles, &Boxes, &Neighbours);

  Store_Parameters (argc, argv, Param);
  printf("stored parameters. Beginning simulation...\n");

  /* main loop of the program */
  State.start_time = clock();
  while (State._time < Param.tf){
    Update_Particles(Param, Particles, State, Neighbours, Boxes);
    State._time += Param.dt;
    Update_Data(Param, Particles, &State, &Data, Neighbours, Boxes);
    Store_Data (Param, Particles, &State, &Data);

    if (State._time > State.NextUpdateProg-EPS){
      Update_Progress(Param, &State);
      State.NextUpdateProg += Param.UpdateInterProg;
    }
  }
  
  // Print run time
  double run_time = (((double)clock()) - ((double)State.start_time)) / CLOCKS_PER_SEC;
  fprintf(Param.output_param,"\nrun_time = %.1f s\n", run_time);
  printf("\nrun_time = %.2f s\n", run_time);

}
