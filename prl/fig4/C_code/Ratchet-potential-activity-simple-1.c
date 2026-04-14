/* 

Similar to Ratchet-potential-activity-simple.c except the following differences:

1. U(x) can have different coefficients b1,b2 from v(x)
2. Patched the issue where terms of v(x) were shifted by the wrong amount.
3. Removed current storage (displacement storage is enough)

*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Library that measures dates and times
#include <assert.h>

#include "./mt19937-64.c" // random number generator

#define EPS 1e-10

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param{
  double dt;                       //Time step
  long N;                          // integers can be int, long, long int.
  double v0;                       // set the scale of the particle speeds
  double a1;                       // coefficient of the 1st fourier mode for activity landscape
  double a2;                       // coefficient of the 2nd fourier mode for activity landscape
  double b1;                       // coefficient of the 1st fourier mode for potential landscape
  double b2;                       // coefficient of the 2nd fourier mode for potential landscape
  int P;                           // number of periods of activity variation that will be included
  double coef1;                    // x coefficient to go inside 1st fourier mode
  double coef2;                    // x coefficient to go inside 2nd fourier mode
  double v0_a1_half;               // v0*a1/2
  double v0_a2_half;               // v0*a2/2
  double epsilon;                  // potential strength
  double potential_coef1;          // coefficient of the potential
  double potential_coef2;          // coefficient of the potential
  double d;                        // phase shift for v
  double shift1;                   // phase shift for v 1st term
  double shift2;                   // phase shift for v 2nd term
  double alpha;                    // tumbling rate
  double alphadt;                  // used to compute the probability of tumbling
  double L;                        // system size
  double final_time;               //end of simulation (real time, not N timesteps)
  double binwidth;                 // width of density/magnetization profile bins
  int Nbin;                        // number of bins for density/magnetization profiles (=L/binwidth)
  int NstepProfile;                // number of steps to use when storing histogram
  double StoreInterProfile;        // Interval between 2 storages of the pos profile
  double StoreInterPos;            // Interval between 2 storages of the position
  double StoreInterDisp;           // Interval between 2 storages of net displacements
} param;
  
// Structure particle contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;     // position of the particle
  double theta; // orientation of particles
  double next_time; // next time to tumble
} particle;

#include "Ratchet-potential-activity-functions-simple-1.c" // our methods


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
  FILE* output_profile_rho;        // File where position profile is stored
  FILE* output_profile_m;          // File where m profile is stored
  FILE* output_pos;                // File where position is stored
  FILE* output_disp;               // File where integrate displacements are stored
  FILE* output_traj;               // File where complete trajectories are stored
  long long seed;                  //  Seed of random number generator

  double NextStoreProfile;         // Next time after which the position profile is stored
  double NextUpdateProfile;        // Next time to add to density/magnetization profile
  double NextStorePos;             // Next time after which the position is stored
  double NextStoreDisp;            // Next time after which the displacement is stored

  
  long* profile_rho;
  long* profile_m;
  
  /* 
     initializing all variables 
  */
  
  Initialize_parameters(argc, argv, 
  &output_param, &output_profile_rho, &output_profile_m, &output_pos, &output_disp, &output_traj, 
  &Param, &_time, &prev_percentage, &Particles, &Displacements, &Integrated_Displacements, &profile_rho, &profile_m, &seed, 
  &NextStoreProfile, &NextUpdateProfile, &NextStorePos, &NextStoreDisp);
  
  init_genrand64(seed);
  
  Store_Parameters(argc, argv, output_param, Param, seed);
  
  start_time = time(NULL);
  
  /* main loop of the program */
  
  while (_time < Param.final_time){
    
    // Move the particles
    Update_Particles(Param,Displacements,Particles,Integrated_Displacements, _time);
    
    // Update time
    _time += Param.dt;
    
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
    }
    
    // If the time is right, store the integrated displacements
    if (_time>NextStoreDisp-EPS){
      Store_Disp(Param,Integrated_Displacements,_time,output_disp);
      NextStoreDisp += Param.StoreInterDisp;
    }
  }
  printf("\n");

  // Free all arrays
  free(Particles);
  free(Displacements);
  free(Integrated_Displacements);
  
  end_time = time(NULL);
  run_time = ((long)end_time) - ((long)start_time);
  printf("run_time = %ld s\n",run_time);
  
}
