#include "./options.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // Library that measures dates and times

#include "../mt19937-64.c" // random number generator

// param is a structure that contains data which are given to
// functions that run the simulations
typedef struct param
{
  double Lx;         // system size in x
  double Ly;         // system size in y
  long N;            // integers can be int, long, long int.
  double dt;         // Time step
  double final_time; // end of simulation (real time, not # steps)

  double T0; // set the scale of the temperature

  // non-Ito discretization
#ifdef ALPHA_NEQ_0
  double alpha;  // 0<alpha<1, discretization timepoint
  double alpha_dt;// alpha*Delta t
#ifdef UNDERDAMP
  raise("Non-Ito discretization not implemented for underdamped particles");
#endif
#ifdef LANDSCAPE_LINEAR
  double Tp_left;   // slope of temperature on left side
  double Tp_right;  // slope of temperature on right side
#endif
#endif

#if !defined(TOFX) & !defined(MUOFX)
#ifdef UNDERDAMP
  double sqrt2gamma0T0dt; // sqrt(2*T0/mu0 *dt)
#elif defined(OVERDAMP)
  double sqrt2D0dt;   // sqrt(2*mu0*T0*dt)
#endif
#endif

#ifdef TOFX
  /*
    Parameters for variation of temperature landscape
  */
  int T_P;          // number of repetitions of motif in x direction
#if defined(LANDSCAPE_CUBIC) || defined(LANDSCAPE_LINEAR)
  double T1;        // temperature at height of variation
  double T_left;    // left boundary of variation
  double T_center1; // left boundary of plateau
  double T_center2; // right boundary of plateau
  double T_right;   // right boundary of variation
  double dT;        // difference between T1 and T0
  double T_dxL;     // spacing between left of variation and left of plateau
  double T_dxR;     // spacing between right of plateau and right of variation
  double T_motif_L; // length of motif
#ifdef LANDSCAPE_CUBIC //cubic interpolation
  double T_a1,T_b1,T_c1,T_d1; // left side interpolation
  double T_a2,T_b2,T_c2,T_d2; // right side interpolation
#endif
#elif defined(LANDSCAPE_TRIG)
  double T_ax1;  // amplitude of 1st mode in x direction
  double T_ax2;  // amplitude of 2nd mode in x direction
  double T_ay1;  // amplitude of 1st mode in y direction
  double T_ay2;  // amplitude of 2nd mode in y direction
  int T_Px;      // number of repetitions of motif in x direction
  int T_Py;      // number of repetitions of motif in y direction
  double T0_ax1; // prefactor of mode
  double T0_ax2; // prefactor of mode
  double T0_ay1; // prefactor of mode
  double T0_ay2; // prefactor of mode
  double T_cx1;  // 2*pi*T_P/L_x
  double T_cx2;  // 4*pi*T_P/L_x
  double T_cy1;  // 2*pi*T_P/L_y
  double T_cy2;  // 4*pi*T_P/L_y
#endif
#endif

  double mu0; // set the scale of the mobility
#ifdef MUOFX
  /*
    Parameters for variation of mobility landscape
    see parameters for temperature landscape for analogous definitions
  */
  int mu_P;       // number of repetitions of motif in x direction
#if defined(LANDSCAPE_CUBIC) || defined(LANDSCAPE_LINEAR)
  double mu1;
  double mu_left;
  double mu_center1;
  double mu_center2;
  double mu_right;
  double dmu;
  double mu_dxL;
  double mu_dxR;
  double mu_motif_L;
#ifdef CUBIC //cubic interpolation
  double mu_a1,mu_b1,mu_c1,mu_d1; // left side interpolation
  double mu_a2,mu_b2,mu_c2,mu_d2; // right side interpolation
#endif
#elif defined(LANDSCAPE_TRIG)
  double mu_ax1;  // amplitude of 1st mode in x direction
  double mu_ax2;  // amplitude of 2nd mode in x direction
  double mu_ay1;  // amplitude of 1st mode in y direction
  double mu_ay2;  // amplitude of 2nd mode in y direction
  double mu_dx;   // phase shift in x direction
  double mu_dy;   // phase shift in y direction
  int mu_Px;      // number of repetitions of motif in x direction
  int mu_Py;      // number of repetitions of motif in y direction
  double mu0_ax1; // prefactor of mode
  double mu0_ax2; // prefactor of mode
  double mu0_ay1; // prefactor of mode
  double mu0_ay2; // prefactor of mode
  double mu_cx1;  // 2*pi*mu_P/L_x
  double mu_cx2;  // 4*pi*mu_P/L_x
  double mu_cy1;  // 2*pi*mu_P/L_y
  double mu_cy2;  // 4*pi*mu_P/L_y
  double mu_cdx1; // 2*pi*dx
  double mu_cdx2; // 4*pi*dx
  double mu_cdy1; // 2*pi*dy
  double mu_cdy2; // 4*pi*dy
#endif
#endif

#if defined(UOFX)
  /*
    Parameters for variation of potential landscape
  */
  int U_P;       // number of repetitions in x direction
#if defined(LANDSCAPE_CUBIC) || defined(LANDSCAPE_LINEAR)
  double U1;        // height of potential hill
  double U_left;    // left boundary of potential hill
  double U_center1; // left boundary of plateau
  double U_center2; // right boundary of plateau
  double U_right;   // right boundary of hill
  double U_dxL;     // spacing between left of variation and left of plateau
  double U_dxR;     // spacing between right of plateau and right of variation
  double U_motif_L; // length of U motif
#ifdef LANDSCAPE_LINEAR
  double F_left;    // force on the left side
  double F_right;   // force on the right side
#elif defined(LANDSCAPE_CUBIC) //cubic interpolation
  double U_a1,U_b1,U_c1; // left side interpolation
  double U_a2,U_b2,U_c2; // right side interpolation
#endif
#elif defined(LANDSCAPE_TRIG)
  double U0;     // overall scale
  double U_ax1;  // amplitude of 1st mode in x direction
  double U_ax2;  // amplitude of 2nd mode in x direction
  double U_ay1;  // amplitude of 1st mode in y direction
  double U_ay2;  // amplitude of 2nd mode in y direction
  double U_dx;   // phase shift in x direction
  double U_dy;   // phase shift in y direction
  int U_Px;      // number of repetitions of motif in x direction
  int U_Py;      // number of repetitions of motif in y direction
  double U0_2pi_ax1; // prefactor of mode
  double U0_4pi_ax2; // prefactor of mode
  double U0_2pi_ay1; // prefactor of mode
  double U0_4pi_ay2; // prefactor of mode
  double U_cx1;  // 2*pi*U_P/L_x
  double U_cx2;  // 4*pi*U_P/L_x
  double U_cy1;  // 2*pi*U_P/L_y
  double U_cy2;  // 4*pi*U_P/L_y
  double U_cdx1; // 2*pi*dx
  double U_cdx2; // 4*pi*dx
  double U_cdy1; // 2*pi*dy
  double U_cdy2; // 4*pi*dy
#endif
#endif

#ifndef NI
  double sigma; // interaction length
  double amp;   // interaction strength
  double rmax2; // interaction cut-off (squared)
  double rmax;  // interaction cutoff
  double rbox;  // width of spatial hashing boxes (spatial units)
  long NxBox;   // number of boxes along x
  long NyBox;   // number of boxes along y
#if defined(FORCE_HARMONIC) || defined(FORCE_INV_HARMONIC)
  double amp_over_sigma;
#elif defined(WCA) || defined(LJ) || defined(FORCE_QUARTIC)
  double force_coef1; // 4*amp*12*sigma^12
  double force_coef2; // 4*amp*6*sigma^6
#endif
#endif

  double binwidthx; // width of position profile bins (only on x axis)
  double binwidthy; // width of position profile bins (only on y axis)
  int Nbinx;        // number of bins for position profiles (=Lx/binwidth)
  int Nbiny;        // number of bins for position profiles (=Ly/binwidth)
#ifdef UNDERDAMP
  double binwidthpx; // width of position profile bins (only on x axis)
  double binwidthpy; // width of position profile bins (only on y axis)
  int Nbinpx;        // number of bins for momentum profiles (=pxmax/binwidth)
  int Nbinpy;        // number of bins for momentum profiles (=pymax/binwidth)
  double pxmax;      // maximum value of momentum along x axis for histogram
  double pymax;      // maximum value of momentum along y axis for histogram
#endif
  double StoreInterProfile; // Interval between 2 storages of the pos
                            // (and, momentum, if underdamped) profile
  int NstepProfile;         // number of steps to use when storing histogram
  double UpdateInterProfile; // Interval between 2 updates of
                             // position(/momentum) profile
  double StoreInterPos;     // Interval between 2 storages of the position
  double StoreInterDisp;    // Interval btwn 2 storages of net position(/momentum)
                            // displacements

#ifdef TRACK
  double Ntrack; // number of particles whose complete trajectory to track
#endif

#ifdef GIVEN_IC
  char *IC_file;
#endif
} param;

// Structure particles which contains all the data needed to characterize the
// state of a particle
typedef struct particle
{
  double x;
  double y; // position of particles
#ifdef EPR
  double x_mid;
  double y_mid; // position of particles at midpoint of timestep
#endif
#ifdef UNDERDAMP
  double px;
  double py; // momentum
#endif
#ifndef NI
  long bi; // box membership
  long bj; // box membership
#ifdef EPR
  long bi_mid; // box membership at midpoint of timestep
  long bj_mid; // box membership at midpoint of timestep
#endif
#endif
} particle;

#ifndef NI
// Structure used to compute forces in the presence of periodic boundary
// conditions
typedef struct box
{
  long i;          // position of box along x axis
  long j;          // position of box along y axis
  double epsilonx; // distance to add along x to take into account PBC
  double epsilony; // distance to add along y to take into account PBC
} box;
#endif

#include "PBPs-2D-functions.c" // our methods

int
main (int argc, char *argv[])
{

  time_t start_time;
  time_t end_time;
  long run_time;

  /* declaration of variables */
  double _time;           // Current simulation time
  double prev_percentage; // previous simulation progress
  particle *Particles;    // Array that contains the state of all particles
  param Param; // Param structure that contains information given to functions
  FILE *output_param;   // File where parameters are stored
  FILE *output_profile; // File where position/momentum profile is stored
  FILE *output_J;       // File where current profile is stored
  FILE *output_pos;     // File where position is stored
  FILE *output_disp;    // File where integrate displacements are stored
#ifdef EPR
  FILE *output_EPR; // File where EPR profile is stored
#endif
#ifdef SIGMA
  FILE* output_profile_sigma;
#endif
#ifdef TRACK
  FILE *output_traj; // File where complete trajectories are stored
#endif
  long long seed; //  Seed of random number generator

  double NextStoreProfile; // Next time after which the density/momentum profile
                           // is stored
  double NextUpdateProfile; // Next time to add to density/momentum profile
  double NextStorePos;      // Next time after which the position is stored
  double NextStoreDisp;     // Next time after which the displacement (in pos. &
                            // mom. space) is stored
#ifdef SIGMA
  int sigma;                // whether or not to update sigma at this point
#endif
  double *Displacements;    // Displacements[2*i] is the displacement along x of
                         // particle i Displacements[2*i+1] is the displacement
                         // along y of particle i.
  double
      *Integrated_Displacements; // Integrated_Displacements[2*i] is the
                                 // cumulative displacement along x on particle
                                 // i Integrated_Displacements[2*i+1] is the
                                 // cumulative displacement along y on particle
                                 // i.
#ifdef UNDERDAMP
  double *Forces; // Forces[2*i] is the force*time increment (p(t+dt) - p(t))
                  // along x of particle i. Forces[2*i+1] is the force*time
                  // increment along y of particle i.
#endif

#ifndef NI
  double *Forces_I; // Forces_I[2*i] is the interaction force along x on
                    // particle i Forces_I[2*i+1] is the interaction force
                    // along y on particle i.
  long **Boxes;     // array containing the first particles in the box.
                    // Box[i][j]=k means the first particle in box i,j is
                    // Particle[k]

  long *Neighbours; // array containing the neighbours of each
                    // particle in a box. Neighbours[2*i+1]=k means
                    // the next particle after i is k.
                    // Neighbours[2*i+1]=-1 means i is the last
                    // particle in its box. Neighbours[2*i]=k
                    // means the particle before i is k.
                    // Neighbours[2*i+1]=-1 means the particle is
                    // the first in the box.

  box ***NeighbouringBoxes; // Neighbouringboxes[i][j] contains the
                            // list of neighbouring boxes of box (i,j)
    
#ifdef EPR          
  double *Forces_I_mid; // interaction forces at middle of timestep (for EPR calculation)
  
  // these are for box/neighbor memberships at MIDPOINT of timestep
  long **Boxes_mid;     // array containing the first particles in the box.
                    // Box[i][j]=k means the first particle in box i,j is
                    // Particle[k]

  long *Neighbours_mid; // array containing the neighbours of each
                    // particle in a box. Neighbours[2*i+1]=k means
                    // the next particle after i is k.
                    // Neighbours[2*i+1]=-1 means i is the last
                    // particle in its box. Neighbours[2*i]=k
                    // means the particle before i is k.
                    // Neighbours[2*i+1]=-1 means the particle is
                    // the first in the box.
#endif
    
#endif

  long ****profile; // density profile
  double ***profile_J;  // current profile
#ifdef EPR
  double **profile_EPR; // EPR profile
#endif
#ifdef SIGMA
  double ***profile_sigma;
#endif

  /*
     initializing all variables
  */
  Initialize_parameters(argc, argv, &output_param, &output_profile, &output_J,
                         &output_pos, &output_disp, 
#ifdef EPR
                         &output_EPR,
#endif
#ifdef SIGMA
                         &output_profile_sigma,
#endif
#ifdef TRACK
                         &output_traj,
#endif
                         &Param, &_time, &prev_percentage, &Particles,
                         &Displacements, &Integrated_Displacements,
#ifdef UNDERDAMP
                         &Forces,
#endif
#ifndef NI
                         &Forces_I, &Boxes, &Neighbours, &NeighbouringBoxes,
#ifdef EPR
                         &Forces_I_mid, &Boxes_mid, &Neighbours_mid,
#endif
#endif
                         &profile, &profile_J,
#ifdef EPR
                         &profile_EPR, 
#endif
#ifdef SIGMA
                         &profile_sigma,
#endif
                         &seed, &NextStoreProfile, 
                         &NextUpdateProfile, 
#ifdef SIGMA
                         &sigma,
#endif
                         &NextStorePos, &NextStoreDisp);

  Store_Parameters (argc, argv, output_param, Param, seed);
  printf("stored parameters. Beginning simulation...\n");

  start_time = time (NULL);

  /* main loop of the program */

  while (_time < Param.final_time)
    {

      // Move the particles
#ifdef NI
#ifdef UNDERDAMP
      Update_Particles_NI (Param, Displacements, Forces, Particles,
                           Integrated_Displacements, 
#ifdef SIGMA
                          sigma, &profile_sigma,
#endif
                          _time);
#else
      Update_Particles_NI (Param, Displacements, Particles,
                           Integrated_Displacements, 
#ifdef SIGMA
                          sigma, &profile_sigma,
#endif
                          _time);
#endif
#else
#ifdef UNDERDAMP
      Update_Particles (Param, Displacements, Forces, Particles,
                        Integrated_Displacements, Forces_I, Neighbours,
                        NeighbouringBoxes, Boxes,
#ifdef SIGMA
                        sigma, &profile_sigma,
#endif
                        _time);
#else
      
      Update_Particles (Param, Displacements, Particles,
                        Integrated_Displacements, Forces_I, Neighbours,
                        NeighbouringBoxes, Boxes,
#ifdef SIGMA
                        sigma, &profile_sigma,
#endif
                        _time);
#endif
#endif
      
#ifdef SIGMA
      // reset sigma (so it won't update unless told to below)
      sigma = 0;
#endif

      // Update time
      _time += Param.dt;

      // For an interval before NextStoreProfile, add to profiles.
      // If NstepProfile=1, we only want to do this when _time=NextStoreProfile,
      // so we ADD epsilon.
      if (_time > NextUpdateProfile - EPS)
        {
          Update_Profile (Param, Particles, &profile);
          Update_Profile_J(Param, Particles, Displacements, forces, &profile_J);
          NextUpdateProfile += Param.UpdateInterProfile;
#ifdef EPR
          Update_EPR (Param, Particles, Displacements, Forces_I_mid, Neighbours_mid, 
                  NeighbouringBoxes, Boxes_mid, &profile_EPR, _time);
#endif
#ifdef SIGMA
          sigma = 1; // turn on update switch (will be passed to Update_Particles(...))
#endif
        }

      // If the time is right, we store the profile
      if (_time > NextStoreProfile - EPS)
        {
          Store_Profile (Param, Particles, &profile, _time, output_profile);
          Store_Profile_J(Param, Particles, &profile_J, _time, output_profile_J);
#ifdef EPR
          Store_EPR (Param, Particles, &profile_EPR, _time, output_EPR);
#endif
#ifdef SIGMA
          Store_Sigma(Param, &profile_sigma, _time, output_profile_sigma);
#endif
          NextStoreProfile += Param.StoreInterProfile;
          NextUpdateProfile
              = NextStoreProfile
                - (Param.NstepProfile - 1) * Param.UpdateInterProfile;

          /*Also Print simulation progresses*/
          PrintSimulationProgress (_time, Param.final_time, &prev_percentage);
        }

      // If the time is right, we store the positions
      if (_time > NextStorePos - EPS)
        {
          Store_Pos (Param, Particles, _time, output_pos);
          NextStorePos += Param.StoreInterPos;

          // If we are storing any complete trajectories, fflush them now.
#ifdef TRACK
          fflush (output_traj);
#endif
        }

      // If the time is right, store the integrated displacements
      if (_time > NextStoreDisp - EPS)
        {
          Store_Disp (Param, Integrated_Displacements, _time, output_disp);
          NextStoreDisp += Param.StoreInterDisp;
        }

#ifdef TRACK
      if (Param.Ntrack > 0)
        {
          Store_Trajectories (Param, Particles, output_traj, _time);
        }
#endif
    }
  printf ("\n");

  // Free all arrays
  free (Particles);
  free (Displacements);
  free (Integrated_Displacements);
#ifdef UNDERDAMP
  free (Forces);
#endif
#ifndef NI
  free (Forces_I);
#endif
  free (profile);

  end_time = time (NULL);
  run_time = ((long)end_time) - ((long)start_time);
  printf ("run_time = %ld s\n", run_time);
}
