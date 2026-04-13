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

#include "../mt19937-64.c" // random number generator

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param{
  double Lx;         // system size in x
  double Ly;         // system size in y
  long N;            // integers can be int, long, long int.
  double dt;         //Time step
  double final_time; //end of simulation (real time, not # steps)

#ifdef RTP
  double alpha;      // tumbling rate
#elif defined ABP
  double Dr;         //rotational diffusivity
  double sqrt2Drdt;  //rotational diffusivity
#endif
  double tau;        // 1/Dr or 1/alpha

#ifdef DT
  double Dt;         // translational diffusivity
  double sqrt2Dtdt;
#endif
  
  double v0;         // set the scale of the particle speeds
  
#ifdef VOFX
  double v1;         // velocity at peak of variation
  double v_center1;  // left side of platform of activity variation (we set v_left=0)
  double v_center2;  // right side of platform of activity variation 
  double v_right;    // right side of the activity variation (we set v_left=0)
  double dv;         // v1 - v0
  double dxR;        // v_right - v_center
  int P;             // number of repititions of activity landscape
  double motif_L;    // length of velocity motif (variation + padding before next one)
#ifdef LANDSCAPE_CUBIC         // cubic interpolation for contibuity
  double c1,d1;      // left interpolation parameters (a1=T0 and b1=0 since v_left=0)
  double a2,b2,c2,d2;// right interpolation parameters
#endif
#endif
  
#ifndef NI
  double sigma;      // interaction length
  double amp;        // interaction strength
  double rmax2;      // interaction cut-off (squared)
  double rmax;       // interaction cutoff
  double rbox;       // width of spatial hashing boxes (spatial units)
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y
#if defined(FORCE_HARMONIC) || defined(FORCE_INV_HARMONIC)
  double amp_over_sigma;
#elif defined(WCA) || defined(LJ) || defined(FORCE_QUARTIC)
  double force_coef1; //4*amp*12*sigma^12 for WCA/LJ or 4*amp/sigma^4 for quartic
  double force_coef2; //4*amp*6*sigma^6 for WCA/LJ or 4*amp*sigma^2 for quartic
#endif
#endif

  double binwidthx;   // width of density/magnetization profile bins (only on x axis)
  double binwidthy;   // width of density/magnetization profile bins (only on y axis)
  int Nbinx;          // number of bins for density/magnetization profiles (=Lx/binwidth)
  int Nbiny;          // number of bins for density/magnetization profiles (=Ly/binwidth)
  int NstepProfile;  // number of steps to use when storing histogram
  double StoreInterProfile;        // Interval between 2 storages of the pos profile
  double StoreInterPos;            // Interval between 2 storages of the position
  double StoreInterDisp;           // Interval between 2 storages of net displacements
#ifdef TRACK
  double Ntrack;     // number of particles whose complete trajectory to track
#endif 

#ifdef GIVEN_IC
  char* IC_file;
#endif
} param;
  

// Structure particles which contains all the data needed to characterize the state of a particle
typedef struct particle{
  double x;
  double y; //position of particles at timestep k
  double theta; // orientation of particles
#ifdef RTP
  double next_time; // next flipping time
#endif
#ifndef NI
  long bi; // position of box along x axis
  long bj; // position of box along y axis
#ifdef EPR
  long bi_mid; // box membership at midpoint of steps k-2 and k-1
  long bj_mid; // box membership at midpoint of steps k-2 and k-1
#endif
#endif
#ifdef EPR
  double x_prev1; // x[k-1]
  double x_prev2; // x[k-2]
  double y_prev1; // y[k-1]
  double y_prev2; // y[k-2]
  double x_mid; // (x[k-1]+x[k-2])/2
  double y_mid; // (y[k-1]+y[k-2])/2
  double dx; // (x[k-1]-x[k-2])
  double dy; 
  double dxdt; // (x[k-1]-x[k-2])/dt
  double dydt; 
  double ddxdt; // (x[k]-2*x[k-1]+x[k-2])/dt (NOT dt^2)
  double ddydt;
#endif
} particle;

#ifndef NI
// Structure used to compute forces in the presence of periodic boundary conditions
typedef struct box{
  long i; // position of box along x axis
  long j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;
#endif

#include "Ratchet-PFAPs-2D-functions.c" // our methods


int main(int argc, char* argv[]){
    
#if defined(EPR) && !defined(ABP)
  printf("EPR only implemented for ABPs.\n");
  return -1;
#endif
#if defined(EPR) && (!defined(FORCE_QUARTIC) && !defined(FORCE_INV_HARMONIC))
  printf("EPR only implemented for quartic and inverse harmonic forces.\n");
  return -1;
#endif
#if defined(NI) && (defined(FORCE_INV_HARMONIC) || defined(HARMONIC) || defined(FORCE_LINEAR) || defined(WCA))
  printf("Must either be interacting or non-interacting!\n");
  return -1;
#endif

  time_t start_time;
  time_t end_time;
  long run_time;

  /* declaration of variables */
  double _time;                    // Current simulation time
  double prev_percentage;          // previous simulation progress
  particle* Particles;             // Array that contains the state of all particles
  param Param;                     // Param structure that contains information given to functions
  FILE* output_param;              // File where parameters are stored
  FILE* output_profile_rho;        // File where position profile is stored
  FILE* output_profile_m;          // File where m profile is stored
  FILE* output_pos;                // File where position is stored
  FILE* output_disp;               // File where integrate displacements are stored
#ifdef EPR
  FILE *output_EPR; // File where EPR profile is stored
#endif
#ifdef SIGMA
  FILE* output_profile_sigma;
  int sigma;               // whether or not to update sigma at this point
#endif
#ifdef TRACK
  FILE* output_traj;               // File where complete trajectories are stored
#endif
#ifdef FORCE_PROFILE
  FILE* output_profile_force;      // File where active force profile is stored
#endif
  long long seed;                  //  Seed of random number generator

  double NextStoreProfile;         // Next time after which the position profile is stored
  double NextUpdateProfile;        // Next time to add to density/magnetization profile
  double NextStorePos;             // Next time after which the position is stored
  double NextStoreDisp;            // Next time after which the displacement is stored
  double* Displacements ;     // Displacements[2*i] is the displacement along x on part. i
	                      // Displacements[2*i+1] is the displacement along y on
	                      // particle i.
  double* Integrated_Displacements ;     // Integrated_Displacements[2*i] is the cumulative displacement along x on part. i
	                                 // Integrated_Displacements[2*i+1] is the cumulative displacement along y on
	                                 // particle i.
  
#ifndef NI
  
  double* forces ;     // forces[2*i] is the force along x on part. i
	               // forces[2*i+1] is the force along y on
	               // particle i.
  long** Boxes;        // array containing the first particles in the
                       // box Box[i][j]=k means the first particle in
                       // box i,j is Particle[k]
  
  long* Neighbours;    // array containing the neighbours of each
		       // particle in a box. Neighbours[2*i+1]=k means
		       // the next particle after i is k. k=-1 means i
		       // is the last particle. Neighbours[2*i]=k
		       // means the particle before i is k. k=-1 means
		       // the particle is the first in the box.

  box *** NeighbouringBoxes; //Neighbouringboxes[i][j] contains the
			     //list of neighbouring boxes of box (i,j)
#endif
  
  // arrays where profiles are saved
  long** profile_rho;
  double** profile_m;
#ifdef EPR
  double ***profile_EPR; // EPR profile
#endif
#ifdef FORCE_PROFILE
  double*** profile_force;
#endif
#ifdef SIGMA
  double*** profile_sigma;
  double*** profile_sigma_dummy; // filler argument for when sigma is not being updated
#endif
  
  /* 
     initializing all variables 
  */
  
  Initialize_parameters(argc, argv, &output_param, &output_profile_rho, &output_profile_m, &output_pos, &output_disp, 
#ifdef TRACK
                         &output_traj, 
#endif
#ifdef FORCE_PROFILE
                         &output_profile_force,
#endif
#ifdef SIGMA
                         &output_profile_sigma,
#endif
                         &Param, &_time, &prev_percentage, &Particles, &Displacements, &Integrated_Displacements, 
#ifndef NI
                         &forces, &Boxes, &Neighbours, &NeighbouringBoxes, 
#endif
                         &profile_rho, &profile_m, 
#ifdef FORCE_PROFILE
                         &profile_force,
#endif
#ifdef SIGMA
                         &profile_sigma,
#endif
                         &seed, &NextStoreProfile, &NextUpdateProfile, &NextStorePos, &NextStoreDisp
#ifdef SIGMA
                         , &sigma
#endif
  );
  
  Store_Parameters(argc, argv, output_param, Param, seed);
  
  start_time = time(NULL);
  
  /* main loop of the program */
  
  while (_time < Param.final_time){
    
#ifdef TRACK
    if (Param.Ntrack>0) {
      Store_Trajectories(Param,Particles, output_traj, _time);
    }
#endif
    
    // Move the particles
#ifdef NI
    Update_Particles_NI(Param,Displacements,Particles,Integrated_Displacements, _time);
#else
    Update_Particles(Param,Displacements,Particles,Integrated_Displacements,forces, 
                      Neighbours, NeighbouringBoxes, Boxes, 
#ifdef SIGMA
                     sigma, &profile_sigma, // sigma may be updating now
#endif
                     _time);
#endif

#ifdef SIGMA
    // reset sigma (so it won't update unless told to below)
    sigma = 0;
#endif
    
    // Update time
    _time += Param.dt;


    
    // For an interval before NextStoreProfile, add to profiles.
    // If NstepProfile=1, we only want to do this when _time=NextStoreProfile, so we ADD epsilon.
    if (_time>NextUpdateProfile-EPS) {
      Update_Profile(Param, Particles, &profile_rho, &profile_m);

#ifdef FORCE_PROFILE
      Update_Force_Profile(Param, Particles, Displacements, &profile_force
#ifndef NI
      , forces
#endif
      );
#endif

#ifdef SIGMA
      sigma = 1; // turn on update switch (will be passed to Update_Particles(...))
#endif
      NextUpdateProfile += Param.tau;
    }



    // If the time is right, we store the profile
    if(_time>NextStoreProfile-EPS){
      Store_Profile(Param, &profile_rho, &profile_m, _time,output_profile_rho,output_profile_m);

#ifdef FORCE_PROFILE
      Store_Force_Profile(Param, &profile_force, _time, output_profile_force);
#endif

#ifdef SIGMA
      Store_Sigma(Param, &profile_sigma, _time, output_profile_sigma);
#endif

      NextStoreProfile += Param.StoreInterProfile;
      NextUpdateProfile = NextStoreProfile - (Param.NstepProfile-1)*Param.tau;
  
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
  }

  
  printf("\n");

  // Free all arrays
  free(Particles);
  free(Displacements);
  free(Integrated_Displacements);
#ifndef NI
  free(forces);
#endif
  free(profile_rho);
  free(profile_m);
  
  end_time = time(NULL);
  run_time = ((long)end_time) - ((long)start_time);
  printf("run_time = %ld s\n",run_time);
  
}
