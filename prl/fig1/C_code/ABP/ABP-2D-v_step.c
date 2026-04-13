// all the imports, and information about the file, here
#include "./options.h"

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param{
  double Lx;                 // system size in x
  double Ly;                 // system size in y
  long N;                    // integers can be int, long, long int.
  double dt;                 // Time step
  double tf;                 // end of simulation (real time, not # steps)
  double tau;                // persistence time

  fparam_lin vParam;         // parameters of activity landscape (included in inhom_fluct.c)
  char* vs;
  char* vxs;

  double sqrt2Drdt;          // for incrementing theta
  
#ifndef NI
  double a;                  // interaction length
  double rmax;               // interaction cutoff
  double k;                  // interaction strength
  double rmax2;              // interaction cut-off (squared)
  long NxBox;                // number of boxes along x
  long NyBox;                // number of boxes along y
  double c1,c2,c3;           // force parameters
  box ***NeighbouringBoxes;
#endif

  double two_pi;             // 2*M_PI
  long long seed;            // Seed of random number generator

  /* Data/parameter storage */


  double UpdateInterProg;    // Interval between 2 updates of progress bar
  double StoreInterPos;      // Interval between 2 storages of the position

  double binwidthx;          // width of density/magnetization profile bins (only on x axis)
  double binwidthy;          // width of density/magnetization profile bins (only on y axis)
  int Nbinx;                 // number of bins for density/magnetization profiles (=Lx/binwidth)
  int Nbiny;                 // number of bins for density/magnetization profiles (=Ly/binwidth)
  double StoreInterProf;     // Interval between 2 storages of the profile
  int NstepProf;             // number of steps to use when storing histogram
  double UpdateInterProf;    // Interval between 2 updates of profile

  FILE* output_param;        // File where parameters are stored
  FILE* output_pos;          // File where position is stored
  FILE* output_prof;         // File where position profile is stored

  double StoreInterDisp;     // Interval between 2 storages of net displacements
  FILE* output_disp;         // File where integrate displacements are stored
} param;
  

// Structure particles which contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;         // x coordinate
  double y;         // y coordinate
  double th;        // orientation of particles
  double dx;        // increment to x
  double dy;        // increment to y
  double ux;        // cos(th)
  double uy;        // sin(th)
  double v_dt;      // self-propulsion speed * dt
  double idx;       // integrated x displacement
  double idy;       // integrated y displacement
#ifndef NI
  double ifx;   // interaction force in x direction
  double ify;   // interaction force in y direction
  int bi;       // box membership
  int bj;       // box membership
#endif
} particle;



// the current state of the simulation
typedef struct state{
  clock_t start_time;              // real time of simulation start
  double _time;                    // Current simulation time
  double prev_percentage;          // previous simulation progress
  double NextUpdateProg;           // next time to print simulation progress

  double NextStorePos;             // Next time after which the position is stored
  double NextStoreProf;            // Next time after which the position profile is stored
  double NextUpdateProf;           // Next time to add to density/magnetization profile
  double NextStoreDisp;            // Next time after which the displacement is stored
} state;


// the accumulated data
typedef struct data{
  long** prof;
} data;

#ifndef NI 
#include "PFs.c"
#endif
#include "./ABP-2D-v_step-functions.c" // our methods

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
  Initialize_Parameters(argc, argv, &Param, &State, &Data, &Particles
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
    Update_Data(Param, Particles, &State, &Data);
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
