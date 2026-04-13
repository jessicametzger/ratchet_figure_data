#include "./options.h"

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param {
  double Lx;                  // system size in x
  double Ly;                  // system size in y
  long N;                     // integers can be int, long, long int.
  double dt;                  // Time step
  double tf;                  // end of simulation (real time, not # steps)
  long long seed;             // Seed of random number generator

#if defined(TOFX_LIN)
  fparam_lin TParam;          // parameters of temperature landscape (included in inhom_fluct.c)
  double sqrt2dt;
  char* Ts;
  char* Txs;
#elif defined(TOFX_CUB)
  fparam_cub TParam;          // parameters of temperature landscape (included in inhom_fluct.c)
  double sqrt2dt;
  char* Ts;
  char* Txs;
#elif defined(TOFX_TRIG)
  char* Tas;
  char* Tbxs;
  char* Tbys;
  fparam_trig TParam;         // parameters of temperature landscape (included in inhom_fluct.c)
  double sqrt2dt;
#else
  double T0;                  // set the scale of the temperature
  double sqrt2T0dt;
#endif

#if defined(FOFX_LIN)
  fparam_lin FParam;          // parameters of force landscape (included in inhom_fluct.c)
  char* Fs;
  char* Fxs;
#elif defined(FOFX_CUB)
  fparam_cub FParam;          // parameters of force landscape (included in inhom_fluct.c)
  char* Fs;
  char* Fxs;
#endif

#ifndef NI
  double a;                   // interaction length
  double rmax;                // interaction cutoff
  double k;                   // interaction strength
  double rmax2;               // interaction cut-off (squared)
  long NxBox;                 // number of boxes along x
  long NyBox;                 // number of boxes along y
  double c1,c2,c3;            // force parameters
  box ***NeighbouringBoxes;
#endif

  double UpdateInterProg;     // Interval between 2 updates of progress bar
  double StoreInterPos;       // Interval between 2 storages of the position
#ifdef STORE_DISP
  double StoreInterDisp;      // Interval btwn 2 storages of net position(/momentum)
                              // displacements
#endif
#if defined(STORE_ANYPROF)
  double binwidthx;           // width of position profile bins (only on x axis)
  double binwidthy;           // width of position profile bins (only on y axis)
  int Nbinx;                  // number of bins for position profiles (=Lx/binwidth)
  int Nbiny;                  // number of bins for position profiles (=Ly/binwidth)
  double StoreInterProf;      // Interval between 2 storages of the pos
                              // (and, momentum, if underdamped) profile
  int NstepProf;              // number of steps to use when storing histogram
  double UpdateInterProf;     // Interval between 2 updates of
                              // position(/momentum) profile
#endif

#ifdef TRACK
  double Ntrack;              // number of particles whose complete trajectory to track
#endif

#ifdef GIVEN_IC
  char *IC_file;
#endif

  FILE *output_param;         // File where parameters are stored
  FILE *output_pos;           // File where position is stored
#ifdef STORE_DISP
  FILE *output_disp;          // File where integrate displacements are stored
#endif
#ifdef STORE_PROF 
  FILE *output_prof;          // File where position/momentum profile is stored
#endif
#ifdef STORE_JPROF
  FILE *output_Jprof;         // File where current profile is stored
#endif
#ifdef STORE_FPROF
  FILE *output_Fprof;         // File where force profile is stored
#endif
#ifdef STORE_TPROF
  FILE *output_Tprof;         // File where temperature profile is stored
#endif
#ifdef STORE_SIGMAIKPROF
  FILE* output_sigmaIKprof;   // file where profile of stress tensor is stored
#endif
#ifdef STORE_EPR
  FILE* output_EPR;
#endif
#ifdef TRACK
  FILE *output_traj;          // File where complete trajectories are stored
#endif
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
#ifdef STORE_JPROF
  double dx_prev; // dx at previous timestep
  double dy_prev; // dy at previous timestep
#endif
#ifdef STORE_DISP
  double idx;   // integrated x displacement
  double idy;   // integrated y displacement
#endif
#ifdef FOFX 
  double Fx;   // x-force * dt
  double Fy;   // y-force * dt
#endif
#ifndef NI
  double ifx;   // interaction force in x direction
  double ify;   // interaction force in y direction
  int bi;       // box membership
  int bj;       // box membership
#ifdef STORE_SIGMAIKPROF
  double* ifxs; // list of x interaction forces due to all others
  double* ifys; // list of y interaction forces due to all others
#endif
#endif
} particle;


typedef struct state{
  clock_t start_time;              // real time of simulation start
  double _time;                    // Current simulation time
  double prev_percentage;          // previous simulation progress
  double NextUpdateProg;           // next time to print simulation progress

  double NextStorePos;             // Next time after which the position is stored
#if defined(STORE_ANYPROF)
  double NextStoreProf;            // Next time after which the position profile is stored
  double NextUpdateProf;           // Next time to add to density/magnetization profile
#endif
#ifdef STORE_DISP
  double NextStoreDisp;            // Next time after which the displacement is stored
#endif
} state;


typedef struct data{
#ifdef STORE_PROF
  long** prof;                     // prof[i][j] = # particles recorded at box x_i,x_j
#endif
#ifdef STORE_SIGMAIKPROF
  double*** sigmaIKprof;           // sigmaIKprof[i][j][k] = stress component k at xi,xj
#endif
#ifdef STORE_JPROF
  double*** Jprof;                 // Jprof[i][j][k] = J component k at xi,xj
#endif
#ifdef STORE_FPROF
  double*** Fprof;                 // Fprof[i][j][k] = rho*F component k at xi,xj
#endif
#ifdef STORE_TPROF
  double** Tprof;                  // Tprof[i][j] = rho*T at xi,xj
#endif
#ifdef STORE_EPR 
  double** EPR;                    // EPR[i][j] = EPR at xi,xj
#endif
} data;

#ifndef NI 
#include "PFs.c"
#endif
#include "./OD-PBPs-2D-functions.c" // our methods

int main (int argc, char *argv[]) {

  /* declaration of variables */
  particle *Particles; // Array that contains the state of all particles
  param Param;         // Param structure that contains information given to functions
  state State;         // state of simulation
  data Data;           // data arrays

#ifndef NI
  long **Boxes;
  long *Neighbours;
#endif

  /* initializing all variables */
  Initialize_parameters(argc, argv, &Param, &State, &Data, &Particles
#ifndef NI
                         , &Boxes, &Neighbours
#endif
                         );

  Store_Parameters (argc, argv, Param);
  printf("stored parameters. Beginning simulation...\n");

  /* main loop of the program */
  State.start_time = clock();
  while (State._time < Param.tf){
    Update_Particles(Param, Particles, State
#ifndef NI
                     , Neighbours, Boxes
#endif
      );
    State._time += Param.dt;
    Update_Data(Param, Particles, &State, &Data
#if !defined(NI) & defined(STORE_EPR)
              , Neighbours, Boxes
#endif
      );
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
