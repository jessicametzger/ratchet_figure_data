/* 
This code simulates N active particles in an activity landscape, to explore the possibility of ratchet.

To measure ratchet currents, we simulate all the particles in parallel, whether they interact or not, and build their statistics until time t_final.

We measure the current by counting the net number of particles that crossed x=0 between measurements.

We store density and magnetization histograms. We occasionally store positions and integrated displacements.
*/


#include <string.h>
#include <stdio.h>
#include "./options.h"
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
#ifdef FORCE_HARMONIC
  double max_amp;    // coefficient of force, = maximum amplitude (=amp/sigma)
#elif defined FORCE_QUARTIC
  double force_coef1;// coefficients of force in quartic model
  double force_coef2;// coefficients of force in quartic model
#elif defined FORCE_LINEAR
  double force;
#endif
  double rmax;       // interaction cutoff
  double rmax2;      // interaction cutoff ^2
  double rbox;       // width of spatial hashing boxes (spatial units)
  int Nbox;          // number of spatial hashing boxes (=L/rbox)
  double binwidth;   // width of density/magnetization profile bins
  int Nbin;          // number of bins for density/magnetization profiles (=L/binwidth)
  int NstepProfile;  // number of steps to use when storing histogram
  double StoreInterCurrent;        // Interval between 2 storages of the current
  double StoreInterProfile;        // Interval between 2 storages of the pos profile
  double StoreInterPos;            // Interval between 2 storages of the position
  double StoreInterDisp;           // Interval between 2 storages of net displacements
#ifdef TRACK
  double Ntrack;     // number of particles whose complete trajectory to track
#elif defined GIVEN_IC
  char* IC_file;
#endif
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
  FILE* output_current;            // File where current is stored
  FILE* output_profile_rho;        // File where position profile is stored
  FILE* output_profile_m;          // File where m profile is stored
  FILE* output_pos;                // File where position is stored
  FILE* output_disp;               // File where integrate displacements are stored
  FILE* output_traj;               // File where complete trajectories are stored
  long long seed;                  //  Seed of random number generator
  
  double* forces;                  // forces[i] is the force along x on particle i

  double NetCrossings;             // Net number of particles that have crossed x=0 in
                                   // current time interval

  double NextStoreCurrent;         // Next time after which the current is stored
  double NextStoreProfile;         // Next time after which the position profile is stored
  double NextUpdateProfile;        // Next time to add to density/magnetization profile
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
  
  long* profile_rho;
  long* profile_m;
  
  /* 
     initializing all variables 
  */
  
  Initialize_parameters(argc, argv, &output_param, &output_current, &output_profile_rho, &output_profile_m, &output_pos, &output_disp, &output_traj, &Param, &_time, &prev_percentage, &Particles, &Displacements, &Integrated_Displacements, &forces, &Boxes, &Neighbors, &NextBoxes, &profile_rho, &profile_m, &seed, &NetCrossings, &NextStoreCurrent, &NextStoreProfile, &NextUpdateProfile, &NextStorePos, &NextStoreDisp);
  
  Store_Parameters(argc, argv, output_param, Param, seed);
  
  start_time = time(NULL);
  
  /* main loop of the program */
  
  while (_time < Param.final_time){
    
    // Move the particles
#ifdef NI
    Update_Particles_NI(Param,Displacements,Particles,Integrated_Displacements, _time, &NetCrossings);
#else
    Update_Particles(Param,Displacements,Particles,Integrated_Displacements,forces, Neighbors, NextBoxes, Boxes, _time, &NetCrossings);
#endif
    
    // Update time
    _time += Param.dt;

    // If the time is right, we store the currents
    if(_time>NextStoreCurrent-EPS){
      Store_Current(Param,NetCrossings,_time,output_current);
      NextStoreCurrent += Param.StoreInterCurrent;
      NetCrossings = 0;
    }
    
    // For an interval before NextStoreProfile, add to profiles.
    // If NstepProfile=1, we only want to do this when _time=NextStoreProfile, so we ADD epsilon.
    if (_time>NextUpdateProfile-EPS) {
      Update_Profile(Param, Particles, profile_rho, profile_m);
      NextUpdateProfile += 2/Param.alpha;
    }

    // If the time is right, we store the profile
    if(_time>NextStoreProfile-EPS){
      Store_Profile(Param,Particles, &profile_rho, &profile_m, _time,output_profile_rho,output_profile_m);
      NextStoreProfile += Param.StoreInterProfile;
      NextUpdateProfile = NextStoreProfile - (Param.NstepProfile-1)*2/Param.alpha;
  
      /*Also Print simulation progresses*/
      PrintSimulationProgress(_time, Param.final_time, &prev_percentage);
    }

    // If the time is right, we store the positions
    if(_time>NextStorePos-EPS){
      Store_Pos(Param,Particles,_time,output_pos);
      NextStorePos += Param.StoreInterPos;
  
      // If we are storing any complete trajectories, fflush them now.
#ifdef TRACK
      fflush(output_traj);
#endif
    }
    
    // If the time is right, store the integrated displacements
    if (_time>NextStoreDisp-EPS){
      Store_Disp(Param,Integrated_Displacements,_time,output_disp);
      NextStoreDisp += Param.StoreInterDisp;
    }
    
#ifdef TRACK
    if (Param.Ntrack>0) {
      Store_Trajectories(Param,Particles, output_traj, _time);
    }
#endif
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
