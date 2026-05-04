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
  
  double v0;         // set the scale of the particle speeds
  
  double v1;         // velocity at peak of variation
  double v_center1;  // left side of platform of activity variation (we set v_left=0)
  double v_center2;  // right side of platform of activity variation 
  double v_right;    // right side of the activity variation (we set v_left=0)
  double dv;         // v1 - v0
  double dxR;        // v_right - v_center
  int P;             // number of repititions of activity landscape
  double motif_L;    // length of velocity motif (variation + padding before next one)
  double c1,d1;      // left interpolation parameters (a1=T0 and b1=0 since v_left=0)
  double a2,b2,c2,d2;// right interpolation parameters
  
#ifndef NI
  double sigma;      // interaction length
  double amp;        // interaction strength
  double rmax2;      // interaction cut-off (squared)
  double rmax;       // interaction cutoff
  double rbox;       // width of spatial hashing boxes (spatial units)
  long NxBox; // number of boxes along x
  long NyBox; // number of boxes along y
  double amp_over_sigma;
#endif

  double binwidthx;   // width of density/magnetization profile bins (only on x axis)
  double binwidthy;   // width of density/magnetization profile bins (only on y axis)
  int Nbinx;          // number of bins for density/magnetization profiles (=Lx/binwidth)
  int Nbiny;          // number of bins for density/magnetization profiles (=Ly/binwidth)
  int NstepProfile;  // number of steps to use when storing histogram
  double StoreInterProfile;        // Interval between 2 storages of the pos profile
  double StoreInterPos;            // Interval between 2 storages of the position
  double StoreInterDisp;           // Interval between 2 storages of net displacements
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
#endif
} particle;

#include "Ratchet-PFAPs-2D-functions.c" // our methods


int main(int argc, char* argv[]){

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
  FILE* output_pos;                // File where position is stored
  FILE* output_disp;               // File where integrate displacements are stored
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
  
  // array where profile is saved
  long** profile_rho;
  
  /* 
     initializing all variables 
  */
  
  Initialize_parameters(argc, argv, &output_param, &output_profile_rho, &output_pos, &output_disp, 
                         &Param, &_time, &prev_percentage, &Particles, &Displacements, &Integrated_Displacements, 
#ifndef NI
                         &forces, &Boxes, &Neighbours, &NeighbouringBoxes, 
#endif
                         &profile_rho,&seed, &NextStoreProfile, &NextUpdateProfile, &NextStorePos, &NextStoreDisp
  );
  
  Store_Parameters(argc, argv, output_param, Param, seed);
  
  start_time = time(NULL);
  
  /* main loop of the program */
  
  while (_time < Param.final_time){
    
    // Move the particles
#ifdef NI
    Update_Particles_NI(Param,Displacements,Particles,Integrated_Displacements, _time);
#else
    Update_Particles(Param,Displacements,Particles,Integrated_Displacements,forces, 
                      Neighbours, NeighbouringBoxes, Boxes, _time);
#endif
    
    // Update time
    _time += Param.dt;


    
    // For an interval before NextStoreProfile, add to profiles.
    // If NstepProfile=1, we only want to do this when _time=NextStoreProfile, so we ADD epsilon.
    if (_time>NextUpdateProfile-EPS) {
      Update_Profile(Param, Particles, &profile_rho);
      NextUpdateProfile += Param.tau;
    }



    // If the time is right, we store the profile
    if(_time>NextStoreProfile-EPS){
      Store_Profile(Param, &profile_rho, _time,output_profile_rho);
      NextStoreProfile += Param.StoreInterProfile;
      NextUpdateProfile = NextStoreProfile - (Param.NstepProfile-1)*Param.tau;
  
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
#ifndef NI
  free(forces);
#endif
  free(profile_rho);
  
  end_time = time(NULL);
  run_time = ((long)end_time) - ((long)start_time);
  printf("run_time = %ld s\n",run_time);
  
}
