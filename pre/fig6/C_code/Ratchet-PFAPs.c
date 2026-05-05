/* 
This code simulates N active particles in an activity landscape, to explore the possibility of ratchet.

To measure ratchet currents, we simulate all the particles in parallel, whether they interact or not, and build their statistics until time t_final.

*/


#define EPS 1e-10

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Library that measures dates and times
#include <assert.h>

#include "./mt19937-64.c" // random number generator

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param{
  double dt;         //Time step
  long N;            // integers can be int, long, long int.
  double v0;         // set the scale of the particle speeds
  double a1;         // coefficient of the 1st fourier mode for activity landscape
  double a2;         // coefficient of the 2nd fourier mode for activity landscape
  int P;             // number of periods of activity variation that will be included
  double coef1;      // x coefficient to go inside 1st fourier mode
  double coef2;      // x coefficient to go inside 2nd fourier mode
  double v0_a1_half; // v0*a1/2
  double v0_a2_half; // v0*a2/2
  double alpha;      // tumbling rate
  double alphadt;    // used to compute the probability of tumbling
  double L;          // system size
  double final_time; //end of simulation (real time, not N timesteps)
  double sigma;      // interaction length
  double amp;        // interaction strength
  double max_amp;    // coefficient of force, = maximum amplitude (=amp/sigma)
  double rmax;       // interaction cutoff
  double rmax2;      // interaction cutoff ^2
  double rbox;       // width of spatial hashing boxes (spatial units)
  int Nbox;          // number of spatial hashing boxes (=L/rbox)
  double StoreInterPos;            // Interval between 2 storages of the position
  double StoreInterDisp;           // Interval between 2 storages of net displacements
} param;
  
// Structure particle contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;     // position of the particle
  double theta; // orientation of particles
  int bi;       // index of its box along x axis
  double next_time;
} particle;

// Used to list a box's neighboring box and shift to account for PBC
typedef struct box{
  int index;        // position of box along x axis
  double epsilonx;  // distance to add along x to take into account PBC
} box;

#include "Ratchet-PFAPs-functions.c" // our methods


int main(int argc, char* argv[]){

  time_t start_time;
  time_t end_time;
  long run_time;

  /* declaration of variables */
  double _time;                    // Current simulation time
  double prev_percentage;          // previous simulation progress
  particle* Particles;             // Array that contains the state of all particles
  param Param;                     // Param structure that contains information given to functions
  double* Displacements;           // Array in which the displacement of all particles is stored
  double* Integrated_Displacements;// Array in which the integrated displacement of particles is stored, from which the current can be computed
  FILE* output_param;              // File where parameters are stored
  FILE* output_pos;                // File where position is stored
  FILE* output_disp;               // File where integrate displacements are stored
  long long seed;                  //  Seed of random number generator
  
  double* forces;                  // forces[i] is the force along x on particle i

  double NextStorePos;             // Next time after which the position is stored
  double NextStoreDisp;            // Next time after which the displacement is stored
  
  long* Boxes;                     // array containing the 1st particle in each box
  long* Neighbors;                 // Neighbors[2*i+1] is the particle after i,
                                   // Neighbors[2*i] is the particle before i.
                                   // If Neighbors[2*i+1]=-1, then i is the last
                                   // particle in the box, and if Neighbors[2*i]=-1, 
                                   // then i is the first particle in the box.
  box* NextBoxes;                  // box after current. NextBoxes[i] is box i+1, or
                                   // box 0 if i=Nbox.
  
  /* 
     initializing all variables 
  */
  
  Initialize_parameters(argc, argv, &output_param, &output_pos, &output_disp, &Param, &_time, &prev_percentage, 
    &Particles, &Displacements, &Integrated_Displacements, &forces, &Boxes, &Neighbors, &NextBoxes, &seed, 
    &NextStorePos, &NextStoreDisp);
  
  Store_Parameters(argc, argv, output_param, Param, seed);
  
  start_time = time(NULL);
  
  /* main loop of the program */
  
  while (_time < Param.final_time){
    
    // Move the particles
    Update_Particles(Param,Displacements,Particles,Integrated_Displacements,forces, Neighbors, NextBoxes, Boxes, _time);
    
    // Update time
    _time += Param.dt;

    // If the time is right, we store the positions
    if(_time>NextStorePos-EPS){
      Store_Pos(Param,Particles,_time,output_pos);
      NextStorePos += Param.StoreInterPos;
    }
    
    // If the time is right, store the integrated displacements
    if (_time>NextStoreDisp-EPS){
      Store_Disp(Param,Integrated_Displacements,_time,output_disp);
      NextStoreDisp += Param.StoreInterDisp;
  
      /*Also Print simulation progresses*/
      PrintSimulationProgress(_time, Param.final_time, &prev_percentage);
    }
  }
  printf("\n");

  // Free all arrays
  free(Particles);
  free(Displacements);
  free(Integrated_Displacements);
  free(forces);
  
  end_time = time(NULL);
  run_time = ((long)end_time) - ((long)start_time);
  printf("run_time = %ld s\n",run_time);
  
}
