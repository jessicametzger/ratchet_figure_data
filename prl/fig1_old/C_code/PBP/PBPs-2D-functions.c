

/*
  To add someone in a box, you need to put it there and to set its new
  neighbours. This functions add at the top.
*/
#ifndef NI
void
AddinBox (long i, long bi, long bj, long **Boxes, long *Neighbours)
{
  long k;
  // Save the old state of Box, -1 or a particle number
  k = Boxes[bi][bj];

  // Put particle i at the top
  Boxes[bi][bj] = i;
  // There is no one above particle i
  Neighbours[2 * i] = -1;
  // k is after i
  Neighbours[2 * i + 1] = k;
  // If k is particle, i is before k
  if (k != -1)
    Neighbours[2 * k] = i;
}

/*
   To remove someone from the box, you have to remove it, but you do
   not need to update its neighbours. This will be done when it is
   added to a new box
*/
void
RemovefromBox (long i, long bi, long bj, long **Boxes, long *Neighbours)
{
  long next;
  long prev;
  // Store the next particle
  next = Neighbours[2 * i + 1];
  // If the particle is the first in the box
  if (Boxes[bi][bj] == i)
    {
      // next is the new top of the box
      Boxes[bi][bj] = next;
      // if next is a particle, there is no one before it
      if (next != -1)
        {
          Neighbours[2 * next] = -1;
        }
    }
  // If there is someone before you
  else
    {
      // Store the previous particle
      prev = Neighbours[2 * i];
      // The next of the previous is your next
      Neighbours[2 * prev + 1] = next;
      // If next is a particle, its previous is your previous
      if (next != -1)
        Neighbours[2 * next] = prev;
    }
}

#ifdef PBC
void
DefineNeighbouringBoxes (box ***NeighbouringBoxes, param Param)
{
  int i, j;
  for (i = 0; i < Param.NxBox; i++)
    {
      NeighbouringBoxes[i] = (box **)malloc (Param.NyBox * sizeof (box *));
      for (j = 0; j < Param.NyBox; j++)
        {
          NeighbouringBoxes[i][j] = (box *)malloc (5 * sizeof (box));

          // This is the same box
          NeighbouringBoxes[i][j][0].i = i;
          NeighbouringBoxes[i][j][0].j = j;
          NeighbouringBoxes[i][j][0].epsilonx = 0;
          NeighbouringBoxes[i][j][0].epsilony = 0;

          // The box above
          NeighbouringBoxes[i][j][1].i = i;
          NeighbouringBoxes[i][j][1].j = (j + 1) % Param.NyBox;
          NeighbouringBoxes[i][j][1].epsilonx = 0;
          NeighbouringBoxes[i][j][1].epsilony = (j == Param.NyBox - 1) ? (Param.Ly) : 0;

          // The box above, to the right
          NeighbouringBoxes[i][j][2].i = (i + 1) % Param.NxBox;
          NeighbouringBoxes[i][j][2].j = (j + 1) % Param.NyBox;
          NeighbouringBoxes[i][j][2].epsilonx = (i == Param.NxBox - 1) ? (Param.Lx) : 0;
          NeighbouringBoxes[i][j][2].epsilony = (j == Param.NyBox - 1) ? (Param.Ly) : 0;

          // The box to the right
          NeighbouringBoxes[i][j][3].i = (i + 1) % Param.NxBox;
          NeighbouringBoxes[i][j][3].j = j;
          NeighbouringBoxes[i][j][3].epsilonx = (i == Param.NxBox - 1) ? (Param.Lx) : 0;
          NeighbouringBoxes[i][j][3].epsilony = 0;

          // The box below, to the right
          NeighbouringBoxes[i][j][4].i = (i + 1) % Param.NxBox;
          NeighbouringBoxes[i][j][4].j = (j - 1 + Param.NyBox) % Param.NyBox;
          NeighbouringBoxes[i][j][4].epsilonx = (i == Param.NxBox - 1) ? (Param.Lx) : 0;
          NeighbouringBoxes[i][j][4].epsilony = (j == 0) ? (-Param.Ly) : 0;
        }
    }
}
#endif

#endif

void
Initialize_parameters (int argc, char *argv[], FILE **output_param, FILE **output_profile,
                       FILE **output_profile_J, FILE **output_pos, FILE **output_disp, 
#ifdef EPR
                       FILE **output_EPR,
#endif
#ifdef SIGMA
                       FILE** output_profile_sigma,
#endif
#ifdef TRACK
                       FILE **output_traj,
#endif
                       param *Param, double *_time, double *prev_percentage, particle **Particles,
                       double **Displacements, double **Integrated_Displacements,
#ifdef UNDERDAMP
                       double **Forces,
#endif
#ifndef NI
                       double **Forces_I, long ***Boxes, long **Neighbours,
                       box ****NeighbouringBoxes,
#ifdef EPR
                       double **Forces_I_mid, long ***Boxes_mid, long **Neighbours_mid,
#endif
#endif
                       long *****profile, double ****profile_J,
#ifdef EPR
                       double ***profile_EPR, 
#endif
#ifdef SIGMA
                       double**** profile_sigma,
#endif
                       long long *seed, double *NextStoreProfile, double *NextUpdateProfile, 
#ifdef SIGMA
                       int* sigma,
#endif
                       double *NextStorePos, double *NextStoreDisp)
{

  long i, j, k, l;              // indices for allocating arrays
  int argctarget = 0;           // Number of parameters that should be used
  char command_base[1000] = ""; // string that contains the desired format of the command line
  char name[200];               // string in which the file names are written
  

  /*
    Format the command line and count the arguments
  */
  argctarget = 0;
  strcat (command_base, "usage: ");
  strcat (command_base, argv[0]); argctarget++;
  strcat (command_base, " ");
  strcat (command_base, "file "); argctarget++;
  strcat (command_base, "Lx "); argctarget++;
  strcat (command_base, "Ly "); argctarget++;
  strcat (command_base, "N "); argctarget++;
  strcat (command_base, "dt "); argctarget++;
  strcat (command_base, "final_time "); argctarget++;
  strcat (command_base, "T0 "); argctarget++;
  strcat (command_base, "mu0 "); argctarget++;
#ifdef ALPHA_NEQ_0
  strcat (command_base, "alpha "); argctarget++;
#endif
#ifdef TOFX
  strcat (command_base, "T_P "); argctarget++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  strcat (command_base, "T1 "); argctarget++;
  strcat (command_base, "T_left "); argctarget++;
  strcat (command_base, "T_center1 "); argctarget++;
  strcat (command_base, "T_center2 "); argctarget++;
  strcat (command_base, "T_right "); argctarget++;
#elif defined(LANDSCAPE_TRIG)
  strcat (command_base, "T_ax1 "); argctarget++;
  strcat (command_base, "T_ax2 "); argctarget++;
  strcat (command_base, "T_ay1 "); argctarget++;
  strcat (command_base, "T_ay2 "); argctarget++;
  strcat (command_base, "T_Py "); argctarget++;
#endif
#endif
    
#ifdef MUOFX
  strcat (command_base, "mu_P "); argctarget++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  strcat (command_base, "mu1 "); argctarget++;
  strcat (command_base, "mu_left "); argctarget++;
  strcat (command_base, "mu_center1 "); argctarget++;
  strcat (command_base, "mu_center2 "); argctarget++;
  strcat (command_base, "mu_right "); argctarget++;
#elif defined(LANDSCAPE_TRIG)
  strcat (command_base, "mu_ax1 "); argctarget++;
  strcat (command_base, "mu_ax2 "); argctarget++;
  strcat (command_base, "mu_ay1 "); argctarget++;
  strcat (command_base, "mu_ay2 "); argctarget++;
  strcat (command_base, "mu_dx "); argctarget++;
  strcat (command_base, "mu_dy "); argctarget++;
  strcat (command_base, "mu_Py "); argctarget++;
#endif
#endif
    
#ifdef UOFX
  strcat (command_base, "U_P "); argctarget++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  strcat (command_base, "U1 "); argctarget++;
  strcat (command_base, "U_left "); argctarget++;
  strcat (command_base, "U_center1 "); argctarget++;
  strcat (command_base, "U_center2 "); argctarget++;
  strcat (command_base, "U_right "); argctarget++;
#elif defined(LANDSCAPE_TRIG)
  strcat (command_base, "U0 "); argctarget++;
  strcat (command_base, "U_ax1 "); argctarget++;
  strcat (command_base, "U_ax2 "); argctarget++;
  strcat (command_base, "U_ay1 "); argctarget++;
  strcat (command_base, "U_ay2 "); argctarget++;
  strcat (command_base, "U_dx "); argctarget++;
  strcat (command_base, "U_dy "); argctarget++;
  strcat (command_base, "U_Py "); argctarget++;
#endif
#endif
    
#ifndef NI
  strcat (command_base, "sigma "); argctarget++;
  strcat (command_base, "amp "); argctarget++;
#endif
  strcat (command_base, "seed "); argctarget++;
  strcat (command_base, "Nbinx "); argctarget++;
  strcat (command_base, "Nbiny "); argctarget++;
#ifdef UNDERDAMP
  strcat (command_base, "Nbinpx "); argctarget++;
  strcat (command_base, "Nbinpy "); argctarget++;
#endif
  strcat (command_base, "StoreInterProfile "); argctarget++;
  strcat (command_base, "NstepProfile "); argctarget++;
  strcat (command_base, "UpdateInterProfile "); argctarget++;
  strcat (command_base, "StoreInterPos "); argctarget++;
  strcat (command_base, "StoreInterDisp "); argctarget++;
#ifdef TRACK
  strcat (command_base, "Ntrack "); argctarget++;
#endif
#ifdef GIVEN_IC
  strcat (command_base, "IC_file "); argctarget++;
#endif

  strcat (command_base, "\n");

  // Check if the call to the program was correct
  if (argc != argctarget)
    {
      printf ("%s\n", command_base);
      exit (1);
    }

  // Affect variables
  i = 1;

  // Read in file name and create files
  sprintf (name, "%s-param", argv[i]);
  output_param[0] = fopen (name, "w");

  sprintf (name, "%s-profile", argv[i]);
  output_profile[0] = fopen (name, "w");

  sprintf (name, "%s-profile_J", argv[i]);
  output_profile_J[0] = fopen (name, "w");

  sprintf (name, "%s-pos", argv[i]);
  output_pos[0] = fopen (name, "w");

  sprintf (name, "%s-disp", argv[i]);
  output_disp[0] = fopen (name, "w");

#ifdef EPR
  sprintf (name, "%s-EPR", argv[i]);
  output_EPR[0] = fopen (name, "w");
#endif
  
#ifdef SIGMA
  sprintf(name,"%s-profile_sigma",argv[i]);
  output_profile_sigma[0]=fopen(name,"w");
#endif

#ifdef TRACK
  sprintf (name, "%s-traj", argv[i]);
  output_traj[0] = fopen (name, "w");
#endif

  printf ("created files\n");

  i++;

  // Read in basic parameters
  Param[0].Lx                 = strtod (argv[i], NULL); i++;
  Param[0].Ly                 = strtod (argv[i], NULL); i++;
  Param[0].N                  = (long)strtod (argv[i], NULL); i++;
  Param[0].dt                 = strtod (argv[i], NULL); i++;
  Param[0].final_time         = strtod (argv[i], NULL); i++;
  Param[0].T0                 = strtod (argv[i], NULL); i++;
  Param[0].mu0                = strtod (argv[i], NULL); i++;
#ifdef ALPHA_NEQ_0
  Param[0].alpha              = strtod (argv[i], NULL); i++;
#endif
#ifdef TOFX
  Param[0].T_P                = (int)strtod (argv[i], NULL); i++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  Param[0].T1                 = strtod (argv[i], NULL); i++;
  Param[0].T_left             = strtod (argv[i], NULL); i++;
  Param[0].T_center1          = strtod (argv[i], NULL); i++;
  Param[0].T_center2          = strtod (argv[i], NULL); i++;
  Param[0].T_right            = strtod (argv[i], NULL); i++;
#elif defined(LANDSCAPE_TRIG)
  Param[0].T_ax1              = strtod (argv[i], NULL); i++;
  Param[0].T_ax2              = strtod (argv[i], NULL); i++;
  Param[0].T_ay1              = strtod (argv[i], NULL); i++;
  Param[0].T_ay2              = strtod (argv[i], NULL); i++;
  Param[0].T_Py               = (int)strtod (argv[i], NULL); i++;
#endif
#endif
#ifdef MUOFX
  Param[0].mu_P               = (int)strtod (argv[i], NULL); i++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  Param[0].mu1                = strtod (argv[i], NULL); i++;
  Param[0].mu_left            = strtod (argv[i], NULL); i++;
  Param[0].mu_center1         = strtod (argv[i], NULL); i++;
  Param[0].mu_center2         = strtod (argv[i], NULL); i++;
  Param[0].mu_right           = strtod (argv[i], NULL); i++;
#elif defined(LANDSCAPE_TRIG)
  Param[0].mu_ax1             = strtod (argv[i], NULL); i++;
  Param[0].mu_ax2             = strtod (argv[i], NULL); i++;
  Param[0].mu_ay1             = strtod (argv[i], NULL); i++;
  Param[0].mu_ay2             = strtod (argv[i], NULL); i++;
  Param[0].mu_dx              = strtod (argv[i], NULL); i++;
  Param[0].mu_dy              = strtod (argv[i], NULL); i++;
  Param[0].mu_Py              = (int)strtod (argv[i], NULL); i++;
#endif
#endif
#ifdef UOFX
  Param[0].U_P                = (int)strtod (argv[i], NULL); i++;
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  Param[0].U1                 = strtod (argv[i], NULL); i++;
  Param[0].U_left             = strtod (argv[i], NULL); i++;
  Param[0].U_center1          = strtod (argv[i], NULL); i++;
  Param[0].U_center2          = strtod (argv[i], NULL); i++;
  Param[0].U_right            = strtod (argv[i], NULL); i++;
#elif defined(LANDSCAPE_TRIG)
  Param[0].U0                 = strtod (argv[i], NULL); i++;
  Param[0].U_ax1              = strtod (argv[i], NULL); i++;
  Param[0].U_ax2              = strtod (argv[i], NULL); i++;
  Param[0].U_ay1              = strtod (argv[i], NULL); i++;
  Param[0].U_ay2              = strtod (argv[i], NULL); i++;
  Param[0].U_dx               = strtod (argv[i], NULL); i++;
  Param[0].U_dy               = strtod (argv[i], NULL); i++;
  Param[0].U_Py               = (int)strtod (argv[i], NULL); i++;
#endif
#endif
#ifndef NI
  Param[0].sigma              = strtod (argv[i], NULL); i++;
  Param[0].amp                = strtod (argv[i], NULL); i++;
#endif
  seed[0]                     = strtod (argv[i], NULL); i++;
  Param[0].Nbinx              = strtod (argv[i], NULL); i++;
  Param[0].Nbiny              = strtod (argv[i], NULL); i++;
#ifdef UNDERDAMP
  Param[0].Nbinpx             = strtod (argv[i], NULL); i++;
  Param[0].Nbinpy             = strtod (argv[i], NULL); i++;
#endif
  Param[0].StoreInterProfile  = strtod (argv[i], NULL); i++;
  Param[0].NstepProfile       = (int)strtod (argv[i], NULL); i++;
  Param[0].UpdateInterProfile = strtod (argv[i], NULL); i++;
  Param[0].StoreInterPos      = strtod (argv[i], NULL); i++;
  Param[0].StoreInterDisp     = strtod (argv[i], NULL); i++;

#ifdef TRACK
  Param[0].Ntrack             = strtod (argv[i], NULL); i++;
#endif
#ifdef GIVEN_IC
  Param[0].IC_file            = argv[i]; i++;
#endif

  // define parameters that are functions of the inputs and/or known

  init_genrand64 (seed[0]);
  _time[0] = 0;
  prev_percentage[0] = 0;

  NextStoreProfile[0] = Param[0].StoreInterProfile;
  NextUpdateProfile[0]
      = NextStoreProfile[0] - (Param[0].UpdateInterProfile) * (Param[0].NstepProfile - 1);
  assert((Param[0].NstepProfile-1)*Param[0].UpdateInterProfile <= Param[0].StoreInterProfile);

#ifdef SIGMA
  sigma[0]                   = 1; // whether or not to update sigma at the next step. start w/yes
#endif
    
  NextStorePos[0] = Param[0].StoreInterPos;
  NextStoreDisp[0] = Param[0].StoreInterDisp;

  assert (Param[0].T0 >= 0.0);
#ifdef TOFX
    
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  assert (Param[0].T1 >= 0.0);
  Param[0].dT = Param[0].T1 - Param[0].T0;
  Param[0].T_motif_L = Param[0].Lx / ((double)Param[0].T_P);

  // landscapes must satisfy left <= center1 <= center2 <= right
  assert (Param[0].T_left <= Param[0].T_center1);
  assert (Param[0].T_center1 <= Param[0].T_center2);
  assert (Param[0].T_center2 <= Param[0].T_right);

  // landscapes must lie inside system
  assert ((Param[0].T_right - Param[0].T_left) * Param[0].T_P <= Param[0].Lx);
  
  Param[0].T_dxL = Param[0].T_center1 - Param[0].T_left;
  Param[0].T_dxR = Param[0].T_right - Param[0].T_center2;
  
#ifdef LANDSCAPE_CUBIC
  // use cubic interpolation from T0 to T1 and vice versa that has derivative zero at the endpoints 
  Param[0].T_a1    = Param[0].dT * pow(Param[0].T_left,2) * (3*Param[0].T_dxL + 2*Param[0].T_left) / pow(Param[0].T_dxL,3) + Param[0].T0;
  Param[0].T_b1    = -6*Param[0].dT * Param[0].T_left * (Param[0].T_dxL + Param[0].T_left) / pow(Param[0].T_dxL,3);
  Param[0].T_c1    = 3*Param[0].dT * (Param[0].T_dxL + 2*Param[0].T_left) / pow(Param[0].T_dxL,3);
  Param[0].T_d1    = -2*Param[0].dT / pow(Param[0].T_dxL,3);
  
  Param[0].T_a2    = -Param[0].dT * pow(Param[0].T_center2,2) * (3*Param[0].T_dxR + 2*Param[0].T_center2) / pow(Param[0].T_dxR,3) + Param[0].T1;
  Param[0].T_b2    = 6*Param[0].dT * Param[0].T_center2 * (Param[0].T_dxR + Param[0].T_center2) / pow(Param[0].T_dxR,3);
  Param[0].T_c2    = -3*Param[0].dT * (Param[0].T_dxR + 2*Param[0].T_center2) / pow(Param[0].T_dxR,3);
  Param[0].T_d2    = 2*Param[0].dT / pow(Param[0].T_dxR,3);
#elif defined(LANDSCAPE_LINEAR) & defined(ALPHA_NEQ_0)
  Param[0].Tp_left = Param[0].dT / Param[0].T_dxL;
  Param[0].Tp_right=-Param[0].dT / Param[0].T_dxR;
  Param[0].alpha_dt = Param[0].dt*Param[0].alpha;
#endif
    
#elif defined(LANDSCAPE_TRIG)
  Param[0].T_Px = Param[0].T_P;
  Param[0].T0_ax1 = Param[0].T0 * Param[0].T_ax1; // prefactor of mode
  Param[0].T0_ax2 = Param[0].T0 * Param[0].T_ax2; // prefactor of mode
  Param[0].T0_ay1 = Param[0].T0 * Param[0].T_ay1; // prefactor of mode
  Param[0].T0_ay2 = Param[0].T0 * Param[0].T_ay2; // prefactor of mode
  Param[0].T_cx1 = 2 * M_PI * Param[0].T_Px / Param[0].Lx;
  Param[0].T_cx2 = 4 * M_PI * Param[0].T_Px / Param[0].Lx;
  Param[0].T_cy1 = 2 * M_PI * Param[0].T_Py / Param[0].Ly;
  Param[0].T_cy2 = 4 * M_PI * Param[0].T_Py / Param[0].Ly;
#endif
    
#elif !defined(MUOFX)
  // constant diffusivity if T,mu constant.
#ifdef UNDERDAMP
  Param[0].sqrt2gamma0T0dt = sqrt (2 * (1 / Param[0].mu0) * Param[0].T0 * dt);
#elif defined(OVERDAMP)
  Param[0].sqrt2D0dt = sqrt (2 * Param[0].T0 * Param[0].mu0 * Param[0].dt);
#endif
#endif

  assert (Param[0].mu0 > 0.0);
#ifdef MUOFX
    
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_LINEAR)
  Param[0].dmu = Param[0].mu1 - Param[0].mu0;
  Param[0].mu_dxL = Param[0].mu_center1 - Param[0].mu_left;
  Param[0].mu_dxR = Param[0].mu_right - Param[0].mu_center2;
  Param[0].mu_motif_L = Param[0].Lx / ((double)Param[0].mu_P);

  // landscapes must satisfy left <= center1 <= center2 <= right
  assert (Param[0].mu_left <= Param[0].mu_center1);
  assert (Param[0].mu_center1 <= Param[0].mu_center2);
  assert (Param[0].mu_center2 <= Param[0].mu_right);

  // landscapes must lie inside system
  assert ((Param[0].mu_right - Param[0].mu_left) * Param[0].mu_P <= Param[0].Lx);
  
#ifdef LANDSCAPE_CUBIC
  // use cubic interpolation from mu0 to mu1 and vice versa that has derivative zero at the endpoints 
  Param[0].mu_a1    = Param[0].dmu * pow(Param[0].mu_left,2) * (3*Param[0].mu_dxL + 2*Param[0].mu_left) / pow(Param[0].mu_dxL,3) + Param[0].mu0;
  Param[0].mu_b1    = -6*Param[0].dmu * Param[0].mu_left * (Param[0].mu_dxL + Param[0].mu_left) / pow(Param[0].mu_dxL,3);
  Param[0].mu_c1    = 3*Param[0].dmu * (Param[0].mu_dxL + 2*Param[0].mu_left) / pow(Param[0].mu_dxL,3);
  Param[0].mu_d1    = -2*Param[0].dmu / pow(Param[0].mu_dxL,3);
  
  Param[0].mu_a2    = -Param[0].dmu * pow(Param[0].mu_center2,2) * (3*Param[0].mu_dxR + 2*Param[0].mu_center2) / pow(Param[0].mu_dxR,3) + Param[0].mu1;
  Param[0].mu_b2    = 6*Param[0].dmu * Param[0].mu_center2 * (Param[0].mu_dxR + Param[0].mu_center2) / pow(Param[0].mu_dxR,3);
  Param[0].mu_c2    = -3*Param[0].dmu * (Param[0].mu_dxR + 2*Param[0].mu_center2) / pow(Param[0].mu_dxR,3);
  Param[0].mu_d2    = 2*Param[0].dmu / pow(Param[0].mu_dxR,3);
#endif
    
#elif defined(LANDSCAPE_TRIG)
  Param[0].mu_Px = Param[0].mu_P;
  Param[0].mu0_ax1 = Param[0].mu0 * Param[0].mu_ax1; // prefactor of mode
  Param[0].mu0_ax2 = Param[0].mu0 * Param[0].mu_ax2; // prefactor of mode
  Param[0].mu0_ay1 = Param[0].mu0 * Param[0].mu_ay1; // prefactor of mode
  Param[0].mu0_ay2 = Param[0].mu0 * Param[0].mu_ay2; // prefactor of mode
  Param[0].mu_cx1 = 2 * M_PI * Param[0].mu_Px / Param[0].Lx;
  Param[0].mu_cx2 = 4 * M_PI * Param[0].mu_Px / Param[0].Lx;
  Param[0].mu_cy1 = 2 * M_PI * Param[0].mu_Py / Param[0].Ly;
  Param[0].mu_cy2 = 4 * M_PI * Param[0].mu_Py / Param[0].Ly;
  Param[0].mu_cdx1 = 2 * M_PI * Param[0].mu_dx;
  Param[0].mu_cdx2 = 4 * M_PI * Param[0].mu_dx;
  Param[0].mu_cdy1 = 2 * M_PI * Param[0].mu_dy;
  Param[0].mu_cdy2 = 4 * M_PI * Param[0].mu_dy;
#endif
#endif

#ifdef UOFX

#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  Param[0].U_dxL = Param[0].U_center1 - Param[0].U_left;
  Param[0].U_dxR = Param[0].U_right - Param[0].U_center2;
  Param[0].U_motif_L = Param[0].Lx / ((double)Param[0].U_P);

  // landscapes must satisfy left <= center1 <= center2 <= right
  assert (Param[0].U_left <= Param[0].U_center1);
  assert (Param[0].U_center1 <= Param[0].U_center2);
  assert (Param[0].U_center2 <= Param[0].U_right);

  // landscapes must lie inside system
  assert ((Param[0].U_right - Param[0].U_left) * Param[0].U_P <= Param[0].Lx);
  
#ifdef LANDSCAPE_LINEAR
  Param[0].F_left = -Param[0].U1 / Param[0].U_dxL; // note: we use U0=0
  Param[0].F_right = Param[0].U1 / Param[0].U_dxR;
#elif defined(LANDSCAPE_CUBIC)
  Param[0].U_a1 = 6*Param[0].U1 * Param[0].U_left * (Param[0].U_dxL + Param[0].U_left) / pow(Param[0].U_dxL,3);
  Param[0].U_b1 = -6*Param[0].U1 * (Param[0].U_dxL + 2*Param[0].U_left) / pow(Param[0].U_dxL,3);
  Param[0].U_c1 = 6*Param[0].U1 / pow(Param[0].U_dxL,3);
   
  Param[0].U_a2 = -6*Param[0].U1 * Param[0].U_center2 * (Param[0].U_dxR + Param[0].U_center2) / pow(Param[0].U_dxR,3);
  Param[0].U_b2 = 6*Param[0].U1 * (Param[0].U_dxR + 2*Param[0].U_center2) / pow(Param[0].U_dxR,3);
  Param[0].U_c2 = -6*Param[0].U1 / pow(Param[0].U_dxR,3);
  
#endif 
    
#elif defined(LANDSCAPE_TRIG)
  Param[0].U_Px = Param[0].U_P;
  Param[0].U0_2pi_ax1 = 2 * M_PI * Param[0].U0 * Param[0].U_ax1 / Param[0].Lx; // prefactor of mode
  Param[0].U0_4pi_ax2 = 4 * M_PI * Param[0].U0 * Param[0].U_ax2 / Param[0].Lx; // prefactor of mode
  Param[0].U0_2pi_ay1 = 2 * M_PI * Param[0].U0 * Param[0].U_ay1 / Param[0].Ly; // prefactor of mode
  Param[0].U0_4pi_ay2 = 4 * M_PI * Param[0].U0 * Param[0].U_ay2 / Param[0].Ly; // prefactor of mode
  Param[0].U_cx1 = 2 * M_PI * Param[0].U_Px / Param[0].Lx;
  Param[0].U_cx2 = 4 * M_PI * Param[0].U_Px / Param[0].Lx;
  Param[0].U_cy1 = 2 * M_PI * Param[0].U_Py / Param[0].Ly;
  Param[0].U_cy2 = 4 * M_PI * Param[0].U_Py / Param[0].Ly;
  Param[0].U_cdx1 = 2 * M_PI * Param[0].U_dx;
  Param[0].U_cdx2 = 4 * M_PI * Param[0].U_dx;
  Param[0].U_cdy1 = 2 * M_PI * Param[0].U_dy;
  Param[0].U_cdy2 = 4 * M_PI * Param[0].U_dy;
#endif
#endif

#ifndef NI
#if defined(FORCE_HARMONIC) || defined(FORCE_INV_HARMONIC)
  Param[0].amp_over_sigma = Param[0].amp / Param[0].sigma;
#elif defined(FORCE_QUARTIC)
  Param[0].force_coef1 = 4*Param[0].amp/pow(Param[0].sigma,4);
  Param[0].force_coef2 = 4*Param[0].amp/pow(Param[0].sigma,2);
#endif
#if defined(WCA)
  Param[0].rmax = Param[0].sigma * pow (2, 1 / 6.);
  Param[0].force_coef1 = 4 * Param.amp * 12 * pow (Param[0].sigma, 12);
  Param[0].force_coef2 = 4 * Param.amp * 6 * pow (Param[0].sigma, 6);
#else
  Param[0].rmax = Param[0].sigma;
#endif
  Param[0].rmax2 = Param[0].rmax * Param[0].rmax;
  Param[0].rbox = Param[0].rmax;
  Param[0].NxBox = (int)(floor (Param[0].Lx / Param[0].rbox) + EPS);
  Param[0].NyBox = (int)(floor (Param[0].Ly / Param[0].rbox) + EPS);

  // System width must be integer multiple of interaction length
  // could change this but this keeps code simple
  assert (abs ((double)Param[0].NxBox - Param[0].Lx / Param[0].rbox) < EPS);
  assert (abs ((double)Param[0].NyBox - Param[0].Ly / Param[0].rbox) < EPS);
#endif

  Param[0].binwidthx = (double)Param[0].Lx / ((double)Param[0].Nbinx);
  Param[0].binwidthy = (double)Param[0].Ly / ((double)Param[0].Nbiny);

#ifdef UNDERDAMP
  // represents the maximum temperature scale
  double T;
#ifdef TOFX
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  T = fmax (Param[0].T0, Param[0].T1);
#elif defined(LANDSCAPE_TRIG)
  T = Param[0].T0 + Param[0].T0_ax1 + Param[0].T0_ax2 + Param[0].T0_bx1 + Param[0].T0_bx2;
#endif
#endif

  // Only keep px, py in histogram from -e*sqrt(T) to e*sqrt(T).
  // This is such that 99% of Boltzmann-distributed points will be
  // accounted for.
  Param[0].pxmax = M_E * sqrt (T);
  Param[0].pymax = M_E * sqrt (T);

  Param[0].binwidthpx = (double)2 * Param[0].pxmax / ((double)Param[0].Nbinpx);
  Param[0].binwidthpy = (double)2 * Param[0].pymax / ((double)Param[0].Nbinpy);
#endif

  // Must have enough time to make enough measurements between profile storages
  assert (Param[0].NstepProfile * Param[0].UpdateInterProfile <= Param[0].StoreInterProfile);

  // Allocate space for simulation data
  Particles[0] = (particle *)malloc (sizeof (particle) * Param[0].N);
  Displacements[0] = (double *)calloc (2 * Param[0].N, sizeof (double));
  Integrated_Displacements[0] = (double *)calloc (2 * Param[0].N, sizeof (double));
#ifdef UNDERDAMP
  Forces[0] = (double *)calloc (2 * Param[0].N, sizeof (double));
#endif
#ifndef NI
  Forces_I[0] = (double *)calloc (2 * Param[0].N, sizeof (double)); 
#ifdef EPR
  Forces_I_mid[0] = (double *)calloc (2 * Param[0].N, sizeof (double)); 
#endif
#endif

  profile[0] = calloc (Param[0].Nbinx, sizeof (long ***));
  for (i = 0; i < Param[0].Nbinx; i++) {
      profile[0][i] = calloc (Param[0].Nbiny, sizeof (long **));
      for (j = 0; j < Param[0].Nbiny; j++) {
#ifdef OVERDAMP
          profile[0][i][j] = calloc (1, sizeof (long *));
          profile[0][i][j][0] = calloc (1, sizeof (long));
#elif defined(UNDERDAMP)
          profile[0][i][j] = calloc (Param[0].Nbinpx, sizeof (long *));
          for (k = 0; k < Param[0].Nbinpx; k++) {
              profile[0][i][j][k] = calloc (Param[0].Nbinpy, sizeof (long));
              for (l=0; l < Param[0].Nbinpy; l++) {
                  profile[0][i][j][k][l] = 0.0;
                }
            }
#endif
        }
    }

  profile_J[0] = calloc(Param[0].Nbinx, sizeof(double**));
  for (i=0; i<Param[0].Nbinx; i++){
    profile_J[0][i] = calloc(Param[0].Nbiny, sizeof(double*));
    for (j=0; j<Param[0].Nbiny; j++){
      profile_J[0][i][j] = calloc(4, sizeof(double));
      profile_J[0][i][j][0] = 0.0; // fx
      profile_J[0][i][j][1] = 0.0; // fy
      profile_J[0][i][j][2] = 0.0; // sqrt(2*D) eta_x
      profile_J[0][i][j][3] = 0.0; // sqrt(2*D) eta_y
    }
  }
    
#ifdef EPR
  profile_EPR[0] = calloc (Param[0].Nbinx, sizeof (long *));
  for (i = 0; i < Param[0].Nbinx; i++)
    {
      profile_EPR[0][i] = calloc (Param[0].Nbiny, sizeof (long ));
      for (j = 0; j < Param[0].Nbiny; j++)
        {
          profile_EPR[0][i][j] = 0.0;
        }
    }
#endif

#ifdef SIGMA
  double Nx, Ny;
#ifdef NI
  Nx = Param[0].Nbinx;
  Ny = Param[0].Nbiny;
#else
  Nx = Param[0].NxBox;
  Ny = Param[0].NyBox;
#endif
  profile_sigma[0]           = (double***) malloc(Nx*sizeof(double**));
  for (i=0;i<Nx;i++){
    profile_sigma[0][i]      = (double**) calloc(Ny,sizeof(double*));
    for (j=0;j<Ny;j++){
#ifdef OVERDAMP
      profile_sigma[0][i][j] = (double*) calloc(4, sizeof(double));
#elif defined(UNDERDAMP)
      profile_sigma[0][i][j] = (double*) calloc(7, sizeof(double));
#endif
      profile_sigma[0][i][j][0] = 0.0; // sigma_xx
      profile_sigma[0][i][j][1] = 0.0; // sigma_xy
      profile_sigma[0][i][j][2] = 0.0; // sigma_yy
      profile_sigma[0][i][j][3] = 0.0; // rho*T
#ifdef UNDERDAMP
      profile_sigma[0][i][j][4] = 0.0; // p_x p_x
      profile_sigma[0][i][j][5] = 0.0; // p_x p_y
      profile_sigma[0][i][j][6] = 0.0; // p_y p_y
#endif
    }
  }
#endif

#ifndef NI
  // Construct list of 1st particle in each box (initialized as empty)
  Boxes[0] = (long **)malloc (Param[0].NxBox * sizeof (long *));
  for (i = 0; i < Param[0].NxBox; i++)
    {
      Boxes[0][i] = (long *)calloc (Param[0].NyBox, sizeof (long));
      for (j = 0; j < Param[0].NyBox; j++)
        Boxes[0][i][j] = -1;
    }
    
  // Construct empty list of neighbors
  Neighbours[0] = (long *)calloc ((Param[0].N + 1) * 2, sizeof (long));
  for (i = 0; i < Param[0].N; i++)
    {
      Neighbours[0][2 * i] = -1;
      Neighbours[0][2 * i + 1] = -1;
    }
    
#ifdef EPR
  Boxes_mid[0] = (long **)malloc (Param[0].NxBox * sizeof (long *));
  for (i = 0; i < Param[0].NxBox; i++)
    {
      Boxes_mid[0][i] = (long *)calloc (Param[0].NyBox, sizeof (long));
      for (j = 0; j < Param[0].NyBox; j++)
        Boxes_mid[0][i][j] = -1;
    }
  Neighbours_mid[0] = (long *)calloc ((Param[0].N + 1) * 2, sizeof (long));
  for (i = 0; i < Param[0].N; i++)
    {
      Neighbours_mid[0][2 * i] = -1;
      Neighbours_mid[0][2 * i + 1] = -1;
    }
#endif

  // Construct list of next boxes
  NeighbouringBoxes[0] = (box ***)malloc (Param[0].NxBox * sizeof (box **));
  DefineNeighbouringBoxes (NeighbouringBoxes[0], Param[0]);
#endif

// initialize particle locations and polarizations
#ifdef RANDOM_IC

#ifdef UNDERDAMP
  double start_angle,
      start_momentum; // initial magnitude/orientation of momentum, randomly
                      // chosen from equilibrium gaussian dist <p^2> = 1/(2T0).
#endif

  for (i = 0; i < Param[0].N; i++)
    {

      // random initialize positions uniformly over interval
      Particles[0][i].x = Param[0].Lx * genrand64_real2 ();
      Particles[0][i].y = Param[0].Ly * genrand64_real2 ();

#ifdef UNDERDAMP
      start_momentum = gasdev () / sqrt (2 * Param[0].T0);
      start_angle = genrand64_real2 ();

      // randomly pick momentum from equilibrium distribution e^(-v^2/(2 T))
      Particles[0][i].px = start_momentum * cos (2 * M_PI * start_angle);
      Particles[0][i].py = start_momentum * sin (2 * M_PI * start_angle);
#endif

#ifndef NI
      // Add particle to its box and link to its neighbors
      Particles[0][i].bi = (int)(floor (Particles[0][i].x / Param[0].rbox) + EPS);
      Particles[0][i].bj = (int)(floor (Particles[0][i].y / Param[0].rbox) + EPS);
      AddinBox (i, Particles[0][i].bi, Particles[0][i].bj, Boxes[0], Neighbours[0]);
#endif
    }

#elif defined GIVEN_IC

  // open file for reading
  FILE *fp;
  fp = fopen (Param[0].IC_file, "r");

  // Read the first line to get the number of columns
  int MAX_LINE_LEN = 1024;
  char line[MAX_LINE_LEN];

  // 4 columns: time, particle ID, x, y, px, py
  int num_columns = 6;
  char *pch;

  // no header
  // Read the last N lines of the file
  // Note: this will throw an error if the number of lines > N
  // It also will only assign the particle ID's that are present in column 2
  // (--> error later if one missing)
  double tmp, x_tmp, y_tmp;
#ifdef UNDERDAMP
  double px_tmp, py_tmp;
#endif
  long id;
  char *ptr;
  while (fgets (line, MAX_LINE_LEN, fp) != NULL)
    {

      // read line to time, id, x, p variables
#ifdef UNDERDAMP
      sscanf (line, "%lg\t%ld\t%lg\t%lg\t%lg\t%lg", &tmp, &id, &x_tmp, &y_tmp, &px_tmp, &py_tmp); 
      Particles[0][id].px = px_tmp;
      Particles[0][id].py = py_tmp;
#else
      sscanf (line, "%lg\t%ld\t%lg\t%lg", &tmp, &id, &x_tmp, &y_tmp);
#endif
      Particles[0][id].x = x_tmp;
      Particles[0][id].y = y_tmp;

#ifndef NI
      // Add particle to its box and link to its neighbors
      Particles[0][id].bi = (int)(floor (Particles[0][id].x / Param[0].rbox) + EPS);
      Particles[0][id].bj = (int)(floor (Particles[0][id].y / Param[0].rbox) + EPS);
      AddinBox (id, Particles[0][id].bi, Particles[0][id].bj, Boxes[0], Neighbours[0]);
#endif
    }

  // set the current time
  _time[0] = tmp;

  NextStoreProfile[0] += _time[0];
  NextUpdateProfile[0] += _time[0];
  NextStorePos[0] += _time[0];
  NextStoreDisp[0] += _time[0];

  // Close the file
  fclose (fp);
#endif

}

#ifdef TOFX
double
Tofx_function (double x,
#ifdef LANDSCAPE_TRIG
               double y,
#endif
               param Param) {
  /*
    Make a piecewise linear temperature field.
    The temperature field will be repeated P times, equally spaced.
  */
    
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in
  n_motif = (int)(floor (x / Param.T_motif_L) + EPS);
  x_inside_motif = x - Param.T_motif_L * ((double)n_motif);
#endif

#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  if (x_inside_motif >= Param.T_right || x_inside_motif <= Param.T_left) return Param.T0;
  else
    {
      if ((x_inside_motif - Param.T_left) < Param.T_center1) {
          return Param.T0 + Param.dT * (x_inside_motif - Param.T_left) / Param.T_dxL;
        }
      else if ((x_inside_motif - Param.T_left) < Param.T_center2) return Param.T1;
      else {
          return Param.T1
                 - Param.dT * ((x_inside_motif - Param.T_left) - Param.T_center2) / Param.T_dxR;
        }
    }
#elif defined(LANDSCAPE_CUBIC)
    // put in piecewise cubic/constant motif as if it is the entire system
  if (x_inside_motif >= Param.T_right || x_inside_motif <= Param.T_left) return Param.T0;
  else
    {
      if ((x_inside_motif - Param.T_left) < Param.T_center1) {
          return Param.T_a1 + Param.T_b1*x_inside_motif + Param.T_c1*pow(x_inside_motif,2) + Param.T_d1*pow(x_inside_motif,3);
        }
      else if ((x_inside_motif - Param.T_left) < Param.T_center2) return Param.T1;
      else {
          return Param.T_a2 + Param.T_b2*x_inside_motif + Param.T_c2*pow(x_inside_motif,2) + Param.T_d2*pow(x_inside_motif,3);
        }
    }
#elif defined(LANDSCAPE_TRIG)
  /*
    Make a trigonometric temperature field.
    The temperature field will be repeated P times, equally spaced.
  */
  return Param.T0 + Param.T0_ax1 * sin (Param.T_cx1 * x)
                  + Param.T0_ax2 * sin (Param.T_cx2 * x)
                  + Param.T0_ay1 * sin (Param.T_cy1 * y) 
                  + Param.T0_ay2 * sin (Param.T_cy2 * y);
#endif
}


#ifdef ALPHA_NEQ_0
double
Tpofx_function (double x,
#ifdef LANDSCAPE_TRIG
               double y,
#endif
               param Param) {
  /*
    derivative of temperature field
  */
    
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in

  n_motif = (int)(floor (x / Param.T_motif_L) + EPS);
  x_inside_motif = x - Param.T_motif_L * ((double)n_motif);
#endif

#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  if (x_inside_motif >= Param.T_right || x_inside_motif <= Param.T_left) return 0.0;
  else
    {
      if ((x_inside_motif - Param.T_left) < Param.T_center1) return Param.Tp_left;
      else if ((x_inside_motif - Param.T_left) < Param.T_center2) return 0.0;
      else return Param.Tp_right;
    }
#elif defined(LANDSCAPE_CUBIC)
  // now, put in piecewise cubic motif as if the motif is the entire system.
  if (x_inside_motif >= Param.T_right || x_inside_motif <= Param.T_left) return 0.0;
  else {
      if ((x_inside_motif - Param.T_left) < Param.T_center1) {
          return Param.T_a1 + Param.T_b1*x_inside_motif + Param.T_c1*pow(x_inside_motif,2);
        }
      else if ((x_inside_motif - Param.T_left) < Param.T_center2) return 0.0;
      else {
          return Param.T_a2 + Param.T_b2*x_inside_motif + Param.T_c2*pow(x_inside_motif,2);
        }
    }
#elif defined(LANDSCAPE_TRIG)
  return Param.T0_2pi_ax1 * cos (Param.T_cx1 * x)
       + Param.T0_4pi_ax2 * cos (Param.T_cx2 * x)
       + Param.T0_2pi_ay1 * cos (Param.T_cy1 * y)
       + Param.T0_4pi_ay2 * cos (Param.T_cy2 * y);
#endif
}
#endif

#endif


#ifdef MUOFX
double
muofx_function (double x, 
#ifdef LANDSCAPE_TRIG
               double y,
#endif
               param Param) {
  /*
    Make a piecewise linear mobility field.
    The mobility field will be repeata1*sin(2*pi*(x + dx1)) + b1*sin(4*pi*(y+dy1))ed P times, equally spaced.
  */
    
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in

  n_motif = (int)(floor (x / Param.T_motif_L) + EPS);
  x_inside_motif = x - Param.T_motif_L * ((double)n_motif);
#endif

#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  if (x_inside_motif >= Param.mu_right || x_inside_motif <= Param.mu_left) return Param.mu0;
  else
    {
      if ((x_inside_motif - Param.mu_left) < Param.mu_center1) {
          return Param.mu0 + Param.dmu * (x_inside_motif - Param.mu_left) / Param.mu_dxL;
        }
      else if ((x_inside_motif - Param.mu_left) < Param.mu_center2) return Param.mu1;
      else {
          return Param.mu1
                 - Param.dmu * ((x_inside_motif - Param.mu_left) - Param.mu_center2) / Param.mu_dxR;
        }
    }
#elif defined(LANDSCAPE_CUBIC)
    // put in piecewise cubic/constant motif as if it is the entire system
  if (x_inside_motif >= Param.mu_right || x_inside_motif <= Param.mu_left) return Param.mu0;
  else
    {
      if ((x_inside_motif - Param.mu_left) < Param.mu_center1) {
          return Param.mu_a1 + Param.mu_b1*x_inside_motif + Param.mu_c1*pow(x_inside_motif,2) + Param.mu_d1*pow(x_inside_motif,3);
        }
      else if ((x_inside_motif - Param.mu_left) < Param.mu_center2) return Param.mu1;
      else {
          return Param.mu_a2 + Param.mu_b2*x_inside_motif + Param.mu_c2*pow(x_inside_motif,2) + Param.mu_d2*pow(x_inside_motif,3);
        }
    }
#elif defined(LANDSCAPE_TRIG)
  return Param.mu0 + Param.mu0_ax1 * sin (Param.mu_cx1 * x + Param.mu_cdx1)
         + Param.mu0_ax2 * sin (Param.mu_cx2 * x + Param.mu_cdx2)
         + Param.mu0_ay1 * sin (Param.mu_cy1 * y + Param.mu_cdy1)
         + Param.mu0_ay2 * sin (Param.mu_cy2 * y + Param.mu_cdy2);
#endif
}
#endif


#ifdef UOFX
double
Fofx_function (double x, 
#ifdef LANDSCAPE_TRIG
               double y,
#endif
               param Param) {

  /*
    Make a piecewise constant force field (piecewise linear potential field).
    The force field will be repeated P times, equally spaced.
  */

#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in

  n_motif = (int)(floor (x / Param.U_motif_L) + EPS);
  x_inside_motif = x - Param.U_motif_L * ((double)n_motif);
#endif

#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  if (x_inside_motif >= Param.U_right || x_inside_motif <= Param.U_left) return 0.0;
  else
    {
      if ((x_inside_motif - Param.U_left) < Param.U_center1) return Param.F_left;
      else if ((x_inside_motif - Param.U_left) < Param.U_center2) return 0.0;
      else return Param.F_right;
    }
#elif defined(LANDSCAPE_CUBIC)
  // now, put in piecewise linear motif as if the motif is the entire system.
  if (x_inside_motif >= Param.U_right || x_inside_motif <= Param.U_left) return 0.0;
  else {
      if ((x_inside_motif - Param.U_left) < Param.U_center1) {
          return Param.U_a1 + Param.U_b1*x_inside_motif + Param.U_c1*pow(x_inside_motif,2);
        }
      else if ((x_inside_motif - Param.U_left) < Param.U_center2) return 0.0;
      else {
          return Param.U_a2 + Param.U_b2*x_inside_motif + Param.U_c2*pow(x_inside_motif,2);
        }
    }
#elif defined(LANDSCAPE_TRIG)
  return Param.U0_2pi_ax1 * cos (Param.U_cx1 * x + Param.U_cdx1)
       + Param.U0_4pi_ax2 * cos (Param.U_cx2 * x + Param.U_cdx2)
       + Param.U0_2pi_ay1 * cos (Param.U_cy1 * y + Param.U_cdy1)
       + Param.U0_4pi_ay2 * cos (Param.U_cy2 * y + Param.U_cdy2);
#endif
}
#endif



#ifndef NI

// Force_harmonic computes the repulsive harmonic force exerted by the particle
// k onto particle j with a concave-down harmonic potential
double
Force_harmonic (double dist, param Param)
{
#ifdef FORCE_HARMONIC
  return Param.amp_over_sigma * dist;
#endif
}

double Force_quartic(double dx, double dx2, param Param){
#ifdef FORCE_QUARTIC
  return -Param.force_coef1*dx2*dx + Param.force_coef2*dx;
#endif
}

// Force_inv_harmonic computes the repulsive harmonic force exerted by the
// particle k onto particle j with a concave-up harmonic potential
double
Force_inv_harmonic (double dist, param Param)
{
#ifdef FORCE_INV_HARMONIC
  return Param.amp - Param.amp_over_sigma * dist;
#endif
}

// compute the repulsive force due to a piecewise linear "tent" interaction
// potential
double
Force_linear (param Param)
{
#ifdef FORCE_LINEAR
  return Param.amp;
#endif
}

// compute the force due to WCA potential (cutoff LJ)
// force_coef1 = 4*amp*12*sigma^12
// force_coef2 = 4*amp*6*sigma^6
void
Force_WCA (double dist, double dist2, param Param)
{
#ifdef WCA
  double dist7 = dist2 * dist2 * dist2 * dist;
  double dist13 = dist7 * dist2 * dist2 * dist2;
  return Param.force_coef1 / dist13 - Param.force_coef2 / dist7;
#endif
}

// This function computes the interactions between particle J and all
// the subsequent particles in the same box
void
Compute_force_same_box (long j, long *Neighbours, double *forces, param Param, particle *Particles, int mid
#ifdef SIGMA
                        ,int sigma, double**** profile_sigma
#endif
                       ) {
  long k;                  // index of neighbors
  double Force_x, Force_y; // Force exerted between two particles
  double dx, dy;           // displacement between particles
  double dist;             // distance between particles
  double dist2;            // distance between particles squared
  double force;            // distance between particles squared
#ifdef SIGMA
  int bi, bj;     // box index
  double sigma_xx, sigma_xy, sigma_yy; // stress tensor components
  bi = Particles[j].bi;
  bj = Particles[j].bj;
#endif

  // Loop through all the following neighbors of j
  k = Neighbours[2 * j + 1]; // start with the 1st neighbor

  // As long as there are neighbors, iterate
  while (k != -1) {
#ifdef EPR
      if (mid){
        dx = Particles[j].x_mid - Particles[k].x_mid; // positive if x_j>x_k
        dy = Particles[j].y_mid - Particles[k].y_mid; // positive if y_j>y_k
      } else {
        dx = Particles[j].x - Particles[k].x; // positive if x_j>x_k
        dy = Particles[j].y - Particles[k].y; // positive if y_j>y_k
      }
#else
      dx = Particles[j].x - Particles[k].x; // positive if x_j>x_k
      dy = Particles[j].y - Particles[k].y; // positive if y_j>y_k
#endif
      dist2 = dx * dx + dy * dy;

      // If the particles are closer than Param.rmax, they interact
      if (dist2 < Param.rmax2)
        {
#ifdef FORCE_HARMONIC
          dist = pow (dist2, 0.5);
          force = Force_harmonic (dist, Param);
#elif defined FORCE_INV_HARMONIC
          dist = pow (dist2, 0.5);
          force = Force_inv_harmonic (dist, Param);
#elif defined FORCE_QUARTIC
          dist = pow (dist2, 0.5);
          force = Force_quartic (dist,dist2,Param);
#elif defined FORCE_LINEAR
          force = Force_linear (Param);
#elif defined WCA
          dist = pow (dist2, 0.5);
          force = Force_WCA (dist, dist2, Param);
#endif
          Force_x = force * dx / dist;
          Force_y = force * dy / dist;
          forces[2 * j] += Force_x;
          forces[2 * j + 1] += Force_y;
          forces[2 * k] -= Force_x;
          forces[2 * k + 1] -= Force_y;
    
#ifdef SIGMA
        // measure stress tensor
        if (sigma){
          sigma_xx = dx*Force_x;
          sigma_xy = dx*Force_y;
          sigma_yy = dy*Force_y;
      
          bi = Particles[j].bi;
          bj = Particles[j].bj;
      
          profile_sigma[0][bi][bj][0] += sigma_xx;
          profile_sigma[0][bi][bj][1] += sigma_xy;
          profile_sigma[0][bi][bj][2] += sigma_yy;
        }
#endif
        }
      k = Neighbours[2 * k + 1]; // k is now the next particle in the list
    }
}

// This functions compute the interactions between particle j and all
// the particles in the box situated on its right, using periodic
// boundary conditions
void
Compute_force_neighbours (long j, int bi, int bj, box ***NeighbouringBoxes, particle *Particles,
                          param Param, double *forces, long *Neighbours, long **Boxes, int mid
#ifdef SIGMA
                          ,int sigma, double**** profile_sigma
#endif
                         ) {
  int nbi, nbj;                   // box name
  int m;                          // neighbour box id
  double dx_box;                  // particle offset
  double dy_box;                  // particle offset
  long k;                         // particle index
  double Force_x, Force_y, force; // Force between particles
  double dx, dy, dist, dist2;     // distance between particles
#ifdef SIGMA
  double sigma_xx, sigma_xy, sigma_yx, sigma_yy; // stress tensor components
  double frac_0, frac_2, frac_4; // frac_k = fraction of separation vector in box k
  double sigma_xx_part, sigma_xy_part, sigma_yy_part; // stress tensor components (piece within box)
  double sigma_xx_part1, sigma_xy_part1, sigma_yy_part1; // stress tensor components (piece within neighbor box)
  double intersect_x; // where it intersects the edge of the box (x coordinate)
  double box_left, box_bottom, box_right, box_top; // coordinates of the box sides
  int nbi1,nbj1; // overlapping box
#endif

  // NeighbouringBoxes[bi][bj][0] is box (bi,bj) itself
  for (m = 1; m < 5; m++) {
      nbi = NeighbouringBoxes[bi][bj][m].i;           // x index of the mth neighbour
      nbj = NeighbouringBoxes[bi][bj][m].j;           // y index of the mth neighbour
      dx_box = NeighbouringBoxes[bi][bj][m].epsilonx; // x offset to be added to the particles in
                                                      // box (bi,bj) to take into account the
                                                      // periodic boundary conditions
      dy_box = NeighbouringBoxes[bi][bj][m].epsilony; // y offset to be added to the particles in
                                                      // box (bi,bj) to take into account the
                                                      // periodic boundary conditions
#ifdef SIGMA
      // get coordinates of edges of box for geometric calculations below
      if (sigma){
        box_left = Param.rbox * (floor(Particles[j].x/Param.rbox)+EPS);
        box_bottom = Param.rbox * (floor(Particles[j].y/Param.rbox)+EPS);
        box_right = box_left + Param.rbox;
        box_top = box_bottom + Param.rbox;
      }
#endif

      // Loop through all particles k in box n
      k = Boxes[nbi][nbj]; // Start with the first

      // As long as there are particles in the box, iterate
      while (k != -1)
        {
#ifdef EPR
          if (mid){
            dx = Particles[j].x_mid - Particles[k].x_mid - dx_box; // positive if x_j>x_k
            dy = Particles[j].y_mid - Particles[k].y_mid - dy_box; // positive if y_j>y_k
          } else {
            dx = Particles[j].x - Particles[k].x - dx_box; // positive if x_j>x_k
            dy = Particles[j].y - Particles[k].y - dy_box; // positive if y_j>y_k
          }
#else
          dx = Particles[j].x - Particles[k].x - dx_box; // positive if x_j>x_k
          dy = Particles[j].y - Particles[k].y - dy_box; // positive if y_j>y_k
#endif
          dist2 = pow (dx, 2) + pow (dy, 2);

          // If the particles are closer than Param.rmax, they interact
          if (dist2 < Param.rmax2) {
#ifdef FORCE_HARMONIC
            dist = pow (dist2, 0.5);
            force = Force_harmonic (dist, Param);
#elif defined FORCE_INV_HARMONIC
            dist = pow (dist2, 0.5);
            force = Force_inv_harmonic (dist, Param);
#elif defined FORCE_QUARTIC
            dist = pow (dist2, 0.5);
            force = Force_quartic (dist,dist2,Param);
#elif defined FORCE_LINEAR
            force = Force_linear (Param);
#elif defined WCA
            dist = pow (dist2, 0.5);
            force = Force_WCA (dist, dist2, Param);
#endif
            Force_x = force * dx / dist;
            Force_y = force * dy / dist;
            forces[2 * j] += Force_x;
            forces[2 * j + 1] += Force_y;
            forces[2 * k] -= Force_x;
            forces[2 * k + 1] -= Force_y;
              
#ifdef SIGMA
            // measure stress tensor
            if (sigma) {

              sigma_xx = dx*Force_x;
              sigma_xy = dx*Force_y;
              sigma_yx = sigma_xy;
              sigma_yy = dy*Force_y;

              if (m==1){ //above
                // dy will be negative
                frac_0 = (box_top-Particles[j].y)/(-dy);

                // add to same box
                sigma_xx_part = sigma_xx*frac_0;
                sigma_xy_part = sigma_xy*frac_0;
                sigma_yy_part = sigma_yy*frac_0;
                
                profile_sigma[0][bi][bj][0] += sigma_xx_part;
                profile_sigma[0][bi][bj][1] += sigma_xy_part;
                profile_sigma[0][bi][bj][2] += sigma_yy_part;

                // add to neighbor box
                profile_sigma[0][nbi][nbj][0] += sigma_xx-sigma_xx_part;
                profile_sigma[0][nbi][nbj][1] += sigma_xy-sigma_xy_part;
                profile_sigma[0][nbi][nbj][2] += sigma_yy-sigma_yy_part;

              } else if (m==2){ //above and to the right
                // dx and dy will be negative
                intersect_x = (box_top - Particles[j].y)*(-dx)/(-dy);
                if (intersect_x > box_right){ // it goes thru box 3
                  frac_0 = (box_right-Particles[j].x)/(-dx);
                  frac_2 = (Particles[k].y+dy_box-box_top)/(-dy);
                  nbi1 = NeighbouringBoxes[bi][bj][3].i; // box 3
                  nbj1 = NeighbouringBoxes[bi][bj][3].j; // box 3
                } else { // goes thru box 1
                  frac_0 = (box_top-Particles[j].y)/(-dy);
                  frac_2 = (Particles[k].x+dx_box-box_right)/(-dx);
                  nbi1 = NeighbouringBoxes[bi][bj][1].i; // box 1
                  nbj1 = NeighbouringBoxes[bi][bj][1].j; // box 1
                }

                // add components (same box)
                sigma_xx_part = sigma_xx*frac_0;
                sigma_xy_part = sigma_xy*frac_0;
                sigma_yy_part = sigma_yy*frac_0;
                
                profile_sigma[0][bi][bj][0] += sigma_xx_part;
                profile_sigma[0][bi][bj][1] += sigma_xy_part;
                profile_sigma[0][bi][bj][2] += sigma_yy_part;

                // add to neighbor box
                sigma_xx_part1 = sigma_xx*frac_2;
                sigma_xy_part1 = sigma_xy*frac_2;
                sigma_yy_part1 = sigma_yy*frac_2;
                profile_sigma[0][nbi][nbj][0] += sigma_xx_part1;
                profile_sigma[0][nbi][nbj][1] += sigma_xy_part1;
                profile_sigma[0][nbi][nbj][2] += sigma_yy_part1;

                // add to overlapping box
                profile_sigma[0][nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
                profile_sigma[0][nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
                profile_sigma[0][nbi1][nbj1][2] += sigma_yy-sigma_yy_part-sigma_yy_part1;

              } else if (m==3){ //to the right
                // dx will be negative
                frac_0 = (box_right-Particles[j].x)/(-dx);

                // add to same box
                sigma_xx_part = sigma_xx*frac_0;
                sigma_xy_part = sigma_xy*frac_0;
                sigma_yy_part = sigma_yy*frac_0;
                
                profile_sigma[0][bi][bj][0] += sigma_xx_part;
                profile_sigma[0][bi][bj][1] += sigma_xy_part;
                profile_sigma[0][bi][bj][2] += sigma_yy_part;

                // add to neighbor box
                profile_sigma[0][nbi][nbj][0] += sigma_xx-sigma_xx_part;
                profile_sigma[0][nbi][nbj][1] += sigma_xy-sigma_xy_part;
                profile_sigma[0][nbi][nbj][2] += sigma_yy-sigma_yy_part;

              } else if (m==4){ //below and to the right
                // dx will be negative and dy will be positive
                intersect_x = (Particles[j].y - box_bottom)*(-dx)/dy;
                if (intersect_x < box_right) { // goes thru box "5"
                  frac_0 = (Particles[j].y - box_bottom)/dy;
                  frac_4 = (Particles[k].x + dx_box - box_right)/(-dx);

                  // will overlap the box below.
                  nbi1 = bi;
                  nbj1 = (bj-1+Param.NyBox)%Param.NyBox;

                } else { // goes thru box 3
                  frac_0 = (box_right - Particles[j].x)/(-dx);
                  frac_4 = (box_bottom - (Particles[k].y+dy_box))/dy;

                  nbi1 = NeighbouringBoxes[bi][bj][3].i;
                  nbj1 = NeighbouringBoxes[bi][bj][3].j;
                }

                // add components (same box)
                sigma_xx_part = sigma_xx*frac_0;
                sigma_xy_part = sigma_xy*frac_0;
                sigma_yy_part = sigma_yy*frac_0;
                
                profile_sigma[0][bi][bj][0] += sigma_xx_part;
                profile_sigma[0][bi][bj][1] += sigma_xy_part;
                profile_sigma[0][bi][bj][2] += sigma_yy_part;

                // add to neighbor box
                sigma_xx_part1 = sigma_xx*frac_4;
                sigma_xy_part1 = sigma_xy*frac_4;
                sigma_yy_part1 = sigma_yy*frac_4;
                profile_sigma[0][nbi][nbj][0] += sigma_xx_part1;
                profile_sigma[0][nbi][nbj][1] += sigma_xy_part1;
                profile_sigma[0][nbi][nbj][2] += sigma_yy_part1;

                // add to overlapping box
                profile_sigma[0][nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
                profile_sigma[0][nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
                profile_sigma[0][nbi1][nbj1][2] += sigma_yy-sigma_yy_part-sigma_yy_part1;
              }
            }
#endif
          }
          k = Neighbours[2 * k + 1]; // k is now the next particle in the list
        }
    }
}

void
Loop_Force (particle *Particles, param Param, double *forces, long *Neighbours,
            box ***NeighbouringBoxes, long **Boxes, int mid
#ifdef SIGMA
            ,int sigma, double**** profile_sigma
#endif
           ) {
  int bi, bj; // Box index
  long j;     // Particle index

  // Initialize the array of forces to zero
  memset (forces, 0, 2 * Param.N * sizeof (double));

  /*
     To compute the force on each particle, we loop through all the boxes and
     compute:
     - the interactions between particles inside the box
     - the interactions between particles inside the box and inside a
     neighboring box
  */

  // Loop through all boxes
  for (bi = 0; bi < Param.NxBox; bi++) {
      for (bj = 0; bj < Param.NyBox; bj++) {

          // Compute the force inside the box

          // Loop through all the particles 'j' in the box.
          j = Boxes[bi][bj]; // Start with j being the 1st particle

          // As long as j is not -1, compute its interactions with all the other
          // particles
          while (j != -1) {

              // Loop first with the particles in the same box
              Compute_force_same_box (j, Neighbours, forces, Param, Particles, mid
#ifdef SIGMA
                                      ,sigma, profile_sigma
#endif
                                     );

              // Loop through the neigboring box
              Compute_force_neighbours (j, bi, bj, NeighbouringBoxes, Particles, Param, forces,
                                        Neighbours, Boxes, mid
#ifdef SIGMA
                                        ,sigma, profile_sigma
#endif
                                       );

              // We have included all the interactions between j and other
              // particles. Now iterate over j
              j = Neighbours[2 * j + 1];
            }
          // We are done with all particles in Box[i], iterate over the box
        }
    }
}

void
Update_Particles (param Param, double *Displacements,
#ifdef UNDERDAMP
                  double *Forces,
#endif
                  particle *Particles, double *Integrated_Displacements,
                  double *Forces_I, long *Neighbours, box ***NeighbouringBoxes, long **Boxes,
#ifdef SIGMA
                  int sigma, double**** profile_sigma,
#endif
                  double _time) {

  long n;
  int bi,bj, newbi, newbj;

  double T;
  double mu;
#ifdef UOFX
  double F;
#endif
#ifdef ALPHA_NEQ_0
  double Tp;
#endif
  double eta_x, eta_y, tmp; // noise

#ifdef UNDERDAMP
  double sqrt2gammaTdt;
#elif defined(OVERDAMP)
  double sqrt2Ddt;
#endif

#if !defined(TOFX) & !defined(MUOFX)
#ifdef UNDERDAMP
  sqrt2gammaTdt = Param.sqrt2gamma0T0dt;
#else
  sqrt2Ddt = Param.sqrt2D0dt;
#endif
#endif

  // incorporate interactions with spatial hashing. 
  // last argument indicates whether or not force is evaluated at midpoint of timestep
  Loop_Force (Particles, Param, Forces_I, Neighbours, NeighbouringBoxes, Boxes, 0
#ifdef SIGMA
              ,sigma, profile_sigma
#endif
             );

  for (n = 0; n < Param.N; n++) {

#ifdef TOFX
      T = Tofx_function(Particles[n].x,
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);

#ifdef ALPHA_NEQ_0
      Tp = Tpofx_function(Particles[n].x,
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);
#endif

#else
      T = Param.T0;
#endif

#ifdef MUOFX
      mu = muofx_function (Particles[n].x, 
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);
#else
      mu = Param.mu0;
#endif

#ifdef UOFX
      F = Fofx_function (Particles[n].x, 
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);
#endif


#if defined(TOFX) || defined(MUOFX)
#ifdef UNDERDAMP
      sqrt2gammaTdt = sqrt (2 * (1 / mu) * T * Param.dt);
#else
      sqrt2Ddt = sqrt (2 * mu * T * Param.dt);
#endif
#endif

      Displacements[2 * n] = 0.0;
      Displacements[2 * n+1] = 0.0;
      
#ifdef UNDERDAMP
      Forces[2 * n] = 0.0;
      Forces[2 * n+1] = 0.0;

      // incorporate interaction forces
      Forces[2 * n] += Forces_I[2 * n] * Param.dt;
      Forces[2 * n + 1] += Forces_I[2 * n + 1] * Param.dt;

#ifdef UOFX
      // external force
      Forces[2 * n] += F * dt;
#endif

      // diffusion
      Forces[2 * n] += sqrt2gammaTdt * gasdev ();
      Forces[2 * n + 1] += sqrt2gammaTdt * gasdev ();

      // drag
      Forces[2*n] += -Param.dt * Particles[n].px / mu;
      Forces[2*n+1] += -Param.dt * Particles[n].py / mu;

      // spatial displacement
      Displacements[2 * n] = Particles[n].px * Param.dt;
      Displacements[2 * n + 1] = Particles[n].py * Param.dt;
//after this, could possibly add other integrators

#else
//overdamped

      eta_x = gasdev();
      eta_y = gasdev();

      // incorporate interactions
      Displacements[2 * n] += Forces_I[2 * n] * Param.dt;
      Displacements[2 * n + 1] += Forces_I[2 * n + 1] * Param.dt;

      // spurious drift
#ifdef ALPHA_NEQ_0
      tmp = Param.alpha_dt*eta_x*Tp;
      Displacements[2*n] += eta_x*tmp;
      Displacements[2*n+1] += eta_y*tmp;
#endif

#ifdef UOFX
      // external force
      Displacements[2 * n] += mu * F * Param.dt;
#endif

      // diffusion
      Displacements[2 * n] += sqrt2Ddt * eta_x;
      Displacements[2 * n + 1] += sqrt2Ddt * eta_x;
#endif


#ifdef SIGMA
      if (sigma){
        bi = Particles[n].bi;
        bj = Particles[n].bj;
        profile_sigma[0][bi][bj][3] += T; // rho_T
#ifdef UNDERDAMP
        profile_sigma[0][bi][bj][4] += Particles[n].px*Particles[n].px; // p_x p_x
        profile_sigma[0][bi][bj][5] += Particles[n].px*Particles[n].py; // p_x p_y
        profile_sigma[0][bi][bj][6] += Particles[n].py*Particles[n].py; // p_y p_y
#endif
      }
#endif
    }

  for (n = 0; n < Param.N; n++) {
      // Once all displacement are computed, move the particles
      Particles[n].x += Displacements[2 * n];
      Particles[n].y += Displacements[2 * n + 1];

#ifdef UNDERDAMP
      Particles[n].px += Forces[2 * n];
      Particles[n].py += Forces[2 * n + 1];
#endif

      // Take care of periodic boundary conditions
#ifdef PBC
      while (Particles[n].x > Param.Lx) Particles[n].x -= Param.Lx;
      while (Particles[n].y > Param.Ly) Particles[n].y -= Param.Ly;
      while (Particles[n].x < 0) Particles[n].x += Param.Lx;
      while (Particles[n].y < 0) Particles[n].y += Param.Ly;
#endif

      // Update box membership
      newbi = (int)(floor (Particles[n].x / Param.rbox) + EPS);
      newbj = (int)(floor (Particles[n].y / Param.rbox) + EPS);

      if (Particles[n].bi != newbi || Particles[n].bj != newbj) {
          RemovefromBox (n, Particles[n].bi, Particles[n].bj, Boxes, Neighbours);
          AddinBox (n, newbi, newbj, Boxes, Neighbours);
          Particles[n].bi = newbi;
          Particles[n].bj = newbj;
        }

      // Update the cumulated displacements
      Integrated_Displacements[2 * n] += Displacements[2 * n];
      Integrated_Displacements[2 * n + 1] += Displacements[2 * n + 1];

    }
}
#endif



#ifdef NI
void
Update_Particles_NI (param Param, double *Displacements,
#ifdef UNDERDAMP
                  double *Forces,
#endif
                  particle *Particles, double *Integrated_Displacements,
#ifdef SIGMA
                  int sigma, double**** profile_sigma,
#endif
                  double _time) {

  long n,j,k;

  double T;
  double mu;
#ifdef UOFX
  double F;
#endif
  double eta_x, eta_y, tmp; // noise

#ifdef UNDERDAMP
  double sqrt2gammaTdt;
#elif defined(OVERDAMP)
  double sqrt2Ddt;
#endif

#if !defined(TOFX) & !defined(MUOFX)
#ifdef UNDERDAMP
  sqrt2gammaTdt = Param.sqrt2gamma0T0dt;
#else
  sqrt2Ddt = Param.sqrt2D0dt;
#endif
#endif

  for (n = 0; n < Param.N; n++) {

#ifdef TOFX
      T = Tofx_function(Particles[n].x,
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);

#ifdef ALPHA_NEQ_0
      Tp = Tpofx_function(Particles[n].x,
#ifdef LANDSCAPE_TRIG
                         Particles[n].y,
#endif
                         Param);
#endif

#else
      T = Param.T0;
#endif

#ifdef MUOFX
      mu = muofx_function (Particles[n].x,
#ifdef LANDSCAPE_TRIG
           Particles[n].y,
#endif
           Param);
#else
      mu = Param.mu0;
#endif

#ifdef UOFX
      F = Fofx_function (Particles[n].x,
#ifdef LANDSCAPE_TRIG
           Particles[n].y,
#endif
           Param);
#endif



#if defined(TOFX) || defined(MUOFX)
#ifdef UNDERDAMP
      sqrt2gammaTdt = sqrt (2 * (1 / mu) * T * Param.dt);
#else
      sqrt2Ddt = sqrt (2 * mu * T * Param.dt);
#endif
#endif

      Displacements[2 * n] = 0.0;
      Displacements[2 * n+1] = 0.0;

#ifdef UNDERDAMP
      Forces[2 * n] = 0.0;
      Forces[2 * n+1] = 0.0;

#ifdef UOFX
      // external force
      Forces[2 * n] += F * dt;
#endif

      // diffusion
      Forces[2 * n] += sqrt2gammaTdt * gasdev ();
      Forces[2 * n + 1] += sqrt2gammaTdt * gasdev ();

      // spatial displacement
      Displacements[2 * n] += Particles[n].px * Param.dt;
      Displacements[2 * n + 1] += Particles[n].py * Param.dt;

#else
      //overdamp

      eta_x = gasdev();
      eta_y = gasdev();

      // spurious drift
#ifdef ALPHA_NEQ_0
      tmp = Param.alpha_dt*eta_x*Tp;
      Displacements[2*n] += eta_x*tmp;
      Displacements[2*n+1] += eta_y*tmp;
#endif

#ifdef UOFX
      // external force
      Displacements[2 * n] += mu * F * Param.dt;
#endif

      // diffusion
      Displacements[2 * n] += sqrt2Ddt * eta_x;
      Displacements[2 * n + 1] += sqrt2Ddt * eta_y;
#endif // OD/UD

#ifdef SIGMA
      if (sigma){
        j = (int)(floor (Particles[n].x / Param.binwidthx) + EPS);
        k = (int)(floor (Particles[n].y / Param.binwidthy) + EPS);
        profile_sigma[0][j][k][3] += T; // rho_T
#ifdef UNDERDAMP
        profile_sigma[0][j][k][4] += Particles[n].px*Particles[n].px; // p_x p_x
        profile_sigma[0][j][k][5] += Particles[n].px*Particles[n].py; // p_x p_y
        profile_sigma[0][j][k][6] += Particles[n].py*Particles[n].py; // p_y p_y
#endif
      }
#endif
      
    }


      for (n = 0; n < Param.N; n++) {

          // Once all displacement are computed, move the particles
          Particles[n].x += Displacements[2 * n];
          Particles[n].y += Displacements[2 * n + 1];

#ifdef UNDERDAMP
          Particles[n].px += Forces[2 * n];
          Particles[n].py += Forces[2 * n + 1];
#endif

          // Take care of periodic boundary conditions
#ifdef PBC
          while (Particles[n].x > Param.Lx) Particles[n].x -= Param.Lx;
          while (Particles[n].y > Param.Ly) Particles[n].y -= Param.Ly;
          while (Particles[n].x < 0) Particles[n].x += Param.Lx;
          while (Particles[n].y < 0) Particles[n].y += Param.Ly;
#endif
          // Update the cumulated displacements
          Integrated_Displacements[2 * n] += Displacements[2 * n];
          Integrated_Displacements[2 * n + 1] += Displacements[2 * n + 1];
        }
}
#endif

  /*
          Print percentage of completed simulation
  */
  void PrintSimulationProgress (double _time, double FinalTime, double *prev_percentage)
  {
    int progress;

    progress = (int)(_time * 100 / FinalTime);
    if (progress >= prev_percentage[0] + 5) {
        printf ("\rIn progress [%d %%]", progress);
        fflush (stdout);
        prev_percentage[0] = progress;
      }
  }

  void Store_Parameters (int argc, char *argv[], FILE *output_param, param Param, long long seed) {
    long i;

    /* Store parameters */
    for (i = 0; i < argc; i++) fprintf (output_param, "%s ", argv[i]);

    fprintf (output_param, "\n");
    fprintf (output_param, "Lx is %lg\n", Param.Lx);
    fprintf (output_param, "Ly is %lg\n", Param.Ly);
    fprintf (output_param, "N is %ld\n", Param.N);
    fprintf (output_param, "dt is %lg\n", Param.dt);
    fprintf (output_param, "final_time is %lg\n", Param.final_time);
    fprintf (output_param, "T0 is %lg\n", Param.T0);
#ifdef ALPHA_NEQ_0
    fprintf(output_param, "alpha is %lg\n", Param.alpha);
#endif
#ifdef TOFX
    fprintf (output_param, "T_P is %d\n", Param.T_P);
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
    fprintf (output_param, "T1 is %lg\n", Param.T1);
    fprintf (output_param, "T_left is %lg\n", Param.T_left);
    fprintf (output_param, "T_center1 is %lg\n", Param.T_center1);
    fprintf (output_param, "T_center2 is %lg\n", Param.T_center2);
    fprintf (output_param, "T_right is %lg\n", Param.T_right);
#elif defined(LANDSCAPE_TRIG)
    fprintf (output_param, "T_ax1 is %lg\n", Param.T_ax1);
    fprintf (output_param, "T_ax2 is %lg\n", Param.T_ax2);
    fprintf (output_param, "T_ay1 is %lg\n", Param.T_ay1);
    fprintf (output_param, "T_ay2 is %lg\n", Param.T_ay2);
    fprintf (output_param, "T_Px is %d\n", Param.T_Px);
    fprintf (output_param, "T_Py is %d\n", Param.T_Py);
#endif
#endif
    fprintf (output_param, "mu0 is %lg\n", Param.mu0);
#ifdef MUOFX
    fprintf (output_param, "mu_P is %d\n", Param.mu_P);
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
    fprintf (output_param, "mu1 is %lg\n", Param.mu1);
    fprintf (output_param, "mu_left is %lg\n", Param.mu_left);
    fprintf (output_param, "mu_center1 is %lg\n", Param.mu_center1);
    fprintf (output_param, "mu_center2 is %lg\n", Param.mu_center2);
    fprintf (output_param, "mu_right is %lg\n", Param.mu_right);
#elif defined(LANDSCAPE_TRIG)
    fprintf (output_param, "mu_ax1 is %lg\n", Param.mu_ax1);
    fprintf (output_param, "mu_ax2 is %lg\n", Param.mu_ax2);
    fprintf (output_param, "mu_ay1 is %lg\n", Param.mu_ay1);
    fprintf (output_param, "mu_ay2 is %lg\n", Param.mu_ay2);
    fprintf (output_param, "mu_dx is %lg\n", Param.mu_dx);
    fprintf (output_param, "mu_dy is %lg\n", Param.mu_dy);
    fprintf (output_param, "mu_Px is %d\n", Param.mu_Px);
    fprintf (output_param, "mu_Py is %d\n", Param.mu_Py);
#endif
#endif
#ifdef UOFX
    fprintf (output_param, "U_P is %d\n", Param.U_P);
#if defined(LANDSCAPE_LINEAR) || defined(LANDSCAPE_CUBIC)
    fprintf (output_param, "U1 is %lg\n", Param.U1);
    fprintf (output_param, "U_left is %lg\n", Param.U_left);
    fprintf (output_param, "U_center1 is %lg\n", Param.U_center1);
    fprintf (output_param, "U_center2 is %lg\n", Param.U_center2);
    fprintf (output_param, "U_right is %lg\n", Param.U_right);
#elif defined(LANDSCAPE_TRIG)
    fprintf (output_param, "U_ax1 is %lg\n", Param.U_ax1);
    fprintf (output_param, "U_ax2 is %lg\n", Param.U_ax2);
    fprintf (output_param, "U_ay1 is %lg\n", Param.U_ay1);
    fprintf (output_param, "U_ay2 is %lg\n", Param.U_ay2);
    fprintf (output_param, "U_dx is %lg\n", Param.U_dx);
    fprintf (output_param, "U_dy is %lg\n", Param.U_dy);
    fprintf (output_param, "U_Px is %d\n", Param.U_Px);
    fprintf (output_param, "U_Py is %d\n", Param.U_Py);
#endif
#endif
#ifndef NI
    fprintf (output_param, "sigma is %lg\n", Param.sigma);
    fprintf (output_param, "amp is %lg\n", Param.amp);
#endif
    fprintf (output_param, "seed is %lld\n", seed);
    fprintf (output_param, "Nbinx is %d\n", Param.Nbinx);
    fprintf (output_param, "Nbiny is %d\n", Param.Nbiny);
#ifdef UNDERDAMP
    fprintf (output_param, "Nbinpx is %d\n", Param.Nbinpx);
    fprintf (output_param, "Nbinpy is %d\n", Param.Nbinpy);
#endif
    fprintf (output_param, "StoreInterProfile is %lg\n", Param.StoreInterProfile);
    fprintf (output_param, "NstepProfile is %d\n", Param.NstepProfile);
    fprintf (output_param, "UpdateInterProfile is %lg\n", Param.UpdateInterProfile);
    fprintf (output_param, "StoreInterPos is %lg\n", Param.StoreInterPos);
    fprintf (output_param, "StoreInterDisp is %lg\n", Param.StoreInterDisp);
#ifdef TRACK
    fprintf (output_param, "Ntrack is %lg\n", Param.Ntrack);
#endif
#ifdef GIVEN_IC
    fprintf (output_param, "IC_file is %lg\n", Param.IC_file);
#endif
    fflush (output_param);
  }

  // Store the identities, positions, and orientations of all particles, with
  // the time
  void Store_Pos (param Param, particle * Particles, double _time, FILE *output_pos) {

    long n;
    for (n = 0; n < Param.N; n++) {
#ifdef UNDERDAMP
        fprintf (output_pos, "%lg\t%ld\t%lg\t%lg\t%lg\t%lg\n", _time, n, Particles[n].x,
                 Particles[n].y, Particles[n].px, Particles[n].py);
#else
      fprintf (output_pos, "%lg\t%ld\t%lg\t%lg\n", _time, n, Particles[n].x, Particles[n].y);
#endif
      }
    fflush (output_pos);
  }

  // Store the identities and integrated displacements of all particles, with
  // the time
  void Store_Disp (param Param, double *Integrated_Displacements, double _time, FILE *output_disp)
  {

    long n;
    for (n = 0; n < Param.N; n++) {
      fprintf (output_disp, "%lg\t%ld\t%lg\t%lg\n", _time, n, Integrated_Displacements[2 * n],
               Integrated_Displacements[2 * n + 1]);
      }
    fflush (output_disp);
  }

  // Update density and magnetization profiles (add one to bin for each particle
  // currently in it)
  void Update_Profile (param Param, particle * Particles, long *****profile)
  {
    long i;   // particle index
    int j, k; // index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
#ifdef UNDERDAMP
    int l, m; // momentum indexto add to. From 0 to Param.Nbinpx-1 and 0 to
              // Param.Nbinpy-1
#endif

    for (i = 0; i < Param.N; i++) {
      j = (int)(floor (Particles[i].x / Param.binwidthx) + EPS);
      k = (int)(floor (Particles[i].y / Param.binwidthy) + EPS);

#ifdef UNDERDAMP

      if (Particles[i].px > Param.pxmax) l = Param.Nbinpx - 1;
      else if (Particles[i].px < -Param.pxmax) l = 0;
      else l = (int)(floor ((Particles[i].px + Param.pxmax) / Param.binwidthpx) + EPS);
      if (Particles[i].py > Param.pymax) m = Param.Nbinpy - 1;
      else if (Particles[i].py < -Param.pymax) m = 0;
      else m = (int)(floor ((Particles[i].py + Param.pymax) / Param.binwidthpy) + EPS);

      profile[0][j][k][l][m] += 1;
#else
      profile[0][j][k][0][0] += 1;
#endif
      }
  }

  void Update_Profile_J(param Param, particle* Particles, double* Displacements, double* forces, double**** profile_J){
    long i;
    int j,k;
    for (i=0; i<Param.N; i++){
      j = (int) (floor(Particles[i].x / Param.binwidthx) + EPS);
      k = (int) (floor(Particles[i].y / Param.binwidthy) + EPS);
      profile_J[0][j][k][0] += forces[2*i];
      profile_J[0][j][k][1] += forces[2*i+1];
      profile_J[0][j][k][2] += Displacements[2*i] - forces[2*i];
      profile_J[0][j][k][3] += Displacements[2*i+1] - forces[2*i+1];
    }
  }

  // Store the profile of positions and polarizations
  void Store_Profile (param Param, particle * Particles, long *****profile, double _time,
                      FILE *output_profile)
  {
    int j, k, l, m; // profile index

    // columns: time, x, y, count
    for (j = 0; j < Param.Nbinx; j++) {
        for (k = 0; k < Param.Nbiny; k++) {
#ifdef UNDERDAMP
            for (l = 0; l < Param.Nbinpx; l++) {
                for (m = 0; m < Param.Nbinpy; m++) {
                    fprintf (output_profile, "%lg\t%lg\t%lg\t%lg\t%lg\t%ld\n", _time,
                             ((double)j) * Param.binwidthx, ((double)k) * Param.binwidthy,
                             -Param.pxmax + ((double)l) * Param.binwidthpx,
                             -Param.pymax + ((double)m) * Param.binwidthpy, profile[0][j][k][l][m]);
                    profile[0][j][k][l][m] = 0;
                  }
              }
#else
          fprintf (output_profile, "%lg\t%lg\t%lg\t%ld\n", _time, ((double)j) * Param.binwidthx,
                   ((double)k) * Param.binwidthy, profile[0][j][k][0][0]);
          profile[0][j][k][0][0] = 0;
#endif
          }
      }
    fflush (output_profile);
  }

void Store_Profile_J(param Param, particle* Particles, double**** profile_J, double _time, FILE* output_profile_J){
  int j,k;
  for (j=0; j<Param.Nbinx; j++){
    for (k=0; k<Param.Nbiny; k++){
      fprintf(output_profile_J, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",_time, ((double)j)*Param.binwidthx,
        ((double)k)*Param.binwidthy, profile_J[0][j][k][0], profile_J[0][j][k][1], profile_J[0][j][k][2], profile_J[0][j][k][3]);
    }
  }
  fflush(output_profile_J);
}
  
#ifdef EPR
  void Update_EPR(param Param, particle * Particles, double *Displacements, double *Forces_I_mid, long *Neighbours_mid, 
                  box ***NeighbouringBoxes, long **Boxes_mid, double ***profile_EPR, double _time){
    long newbi_mid, newbj_mid;
    int j,k; // profile index
    long n; // particle index
    double epr_tmp, T_tmp, mu_tmp;
#ifdef SIGMA
    double ***dummy;
#endif

    // calculate midpoint data
    for (n=0; n<Param.N; n++){
        
      // move back my half a timestep
      Particles[n].x_mid = Particles[n].x - Displacements[2*n]/2.0;
      Particles[n].y_mid = Particles[n].y - Displacements[2*n+1]/2.0;
    
      // Take care of periodic boundary conditions
#ifdef PBC
      while (Particles[n].x_mid > Param.Lx) Particles[n].x_mid -= Param.Lx;
      while (Particles[n].y_mid > Param.Ly) Particles[n].y_mid -= Param.Ly;
      while (Particles[n].x_mid < 0) Particles[n].x_mid += Param.Lx;
      while (Particles[n].y_mid < 0) Particles[n].y_mid += Param.Ly;
#endif
    
#ifndef NI
      // Update box membership
      newbi_mid = (int)(floor (Particles[n].x_mid / Param.rbox) + EPS);
      newbj_mid = (int)(floor (Particles[n].y_mid / Param.rbox) + EPS);

      if (Particles[n].bi_mid != newbi_mid || Particles[n].bj_mid != newbj_mid)
        {
          RemovefromBox (n, Particles[n].bi_mid, Particles[n].bj_mid, Boxes_mid, Neighbours_mid);
          AddinBox (n, newbi_mid, newbj_mid, Boxes_mid, Neighbours_mid);
          Particles[n].bi_mid = newbi_mid;
          Particles[n].bj_mid = newbj_mid;
        }
      }
    
    // incorporate interactions with spatial hashing
    // last argument indicates whether or not force is evaluated at midpoint of timestep
    Loop_Force (Particles, Param, Forces_I_mid, Neighbours_mid, NeighbouringBoxes, Boxes_mid, 1
#ifdef SIGMA
                ,0, &dummy // "dummy" is here just to fill the argument. It will not be updated since sigma=0.
#endif
               );
#endif
    

    for (n=0; n<Param.N; n++){
      
      j = (int)(floor (Particles[n].x_mid / Param.binwidthx) + EPS);
      k = (int)(floor (Particles[n].y_mid / Param.binwidthy) + EPS);
      
      epr_tmp = 0.0;
#ifdef TOFX
      T_tmp = Tofx_function(Particles[n].x_mid, 
#ifdef LANDSCAPE_TRIG
                            Particles[n].y_mid,
#endif
                            Param);
#else
      T_tmp = Param.T0;
#endif
      
#ifdef MUOFX
      mu_tmp = muofx_function(Particles[n].x_mid, 
#ifdef LANDSCAPE_TRIG
                            Particles[n].y_mid,
#endif
                            Param);
#else
      mu_tmp = Param.mu0;
#endif
      
#if defined(UOFX)
      epr_tmp += Fofx_function(Particles[n].x_mid,
#ifdef LANDSCAPE_TRIG
           Particles[n].y,
#endif
           Param) * Displacements[2*n] / (T_tmp/mu_tmp);
#endif
      
#ifndef NI
      epr_tmp += (Forces_I_mid[2*n] * Displacements[2*n] + Forces_I_mid[2*n+1] * Displacements[2*n+1]) / (T_tmp/mu_tmp);
#endif
#ifdef UNDERDAMP
      epr_tmp -= (Forces[2*n] * Displacements[2*n] + Forces[2*n+1] * Displacements[2*n+1]) / (Param.dt*T_tmp);
#endif
      
      profile_EPR[0][j][k] += epr_tmp;
    
      }
  }
  
  void Store_EPR(param Param, particle *Particles, double ***profile_EPR, double _time, FILE *output_EPR){
    int j,k;
    for (j = 0; j < Param.Nbinx; j++) {
        for (k = 0; k < Param.Nbiny; k++) {
          //columns: time, x, y, epr
          fprintf(output_EPR, "%lg\t%lg\t%lg\t%lg\n", _time, ((double)j)*Param.binwidthx,
                  ((double)k)*Param.binwidthy, profile_EPR[0][j][k]);
          profile_EPR[0][j][k]=0.0;
          }
      }
  }
#endif


#ifdef SIGMA
void Store_Sigma(param Param, double**** profile_sigma, double _time, FILE* output_profile_sigma) {
  int j,k;
  int Nx,Ny; // number of x, y boxes
  double lxbox,lybox; // width of x,y boxes
#ifdef NI
  Nx = Param.Nbinx;
  Ny = Param.Nbiny;
  lxbox = Param.binwidthx;
  lybox = Param.binwidthy;
#else
  Nx = Param.NxBox;
  Ny = Param.NyBox;
  lxbox = Param.rbox;
  lybox = Param.rbox;
#endif
  
  for (j=0; j<Nx; j++) {
    for (k=0; k<Ny; k++) {

#ifdef UNDERDAMP
      fprintf(output_profile_sigma, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", _time, 
             ((double)j)*lxbox, ((double)k)*lybox,
             profile_sigma[0][j][k][0], profile_sigma[0][j][k][1], 
             profile_sigma[0][j][k][2], profile_sigma[0][j][k][3],
             profile_sigma[0][j][k][4], profile_sigma[0][j][k][5], 
             profile_sigma[0][j][k][6]);
      profile_sigma[0][j][k][4]=0.0; // p_x p_x
      profile_sigma[0][j][k][5]=0.0; // p_x p_y
      profile_sigma[0][j][k][6]=0.0; // p_y p_y
#else
      fprintf(output_profile_sigma, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", _time, 
             ((double)j)*lxbox, ((double)k)*lybox,
             profile_sigma[0][j][k][0], profile_sigma[0][j][k][1], 
             profile_sigma[0][j][k][2], profile_sigma[0][j][k][3]);
#endif
      profile_sigma[0][j][k][0]=0.0; // sigmaik_xx
      profile_sigma[0][j][k][1]=0.0; // sigmaik_xy
      profile_sigma[0][j][k][2]=0.0; // sigmaik_yy
      profile_sigma[0][j][k][3]=0.0; // rho*T
    }
  }
  fflush(output_profile_sigma);
}
#endif



#ifdef TRACK
  // Store the complete trajectories (position and polarization) of first Ntrack
  // particles
  void Store_Trajectories (param Param, particle * Particles, FILE * output_traj, double _time)
  {
    long n; // particle index

    for (n = 0; n < Param.Ntrack; n++) {
#ifdef UNDERDAMP
        fprintf (output_traj, "%lg\t%ld\t%lg\t%lg\t%lg\t%lg\t%lg\n", _time, n, Particles[n].x,
                 Particles[n].y, Particles[n].px, Particles[n].py);
#else
        fprintf (output_traj, "%lg\t%ld\t%lg\t%lg\n", _time, n, Particles[n].x,
                 Particles[n].y);
#endif
      }
  }
#endif
