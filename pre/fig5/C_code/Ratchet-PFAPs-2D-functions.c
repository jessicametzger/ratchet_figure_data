

#ifndef NI
#include "/home/jessica/Code/util/spatial_hashing.c"
#endif


void Initialize_parameters(int argc, char* argv[],
  FILE** output_param, FILE** output_profile_rho, FILE** output_profile_m, FILE** output_pos, FILE** output_disp, 
#ifdef TRACK
  FILE** output_traj, 
#endif
#ifdef FORCE_PROFILE
  FILE** output_profile_force,
#endif
#ifdef SIGMA
  FILE** output_profile_sigma,
#endif
  param* Param, double* _time, double* prev_percentage,
  particle** Particles,double** Displacements, double** Integrated_Displacements, 
#ifndef NI
  double** forces, long*** Boxes, long** Neighbours, box**** NeighbouringBoxes,
#endif
  long*** profile_rho, double*** profile_m,
#ifdef FORCE_PROFILE
  double**** profile_force,
#endif
#ifdef SIGMA
  double**** profile_sigma,
#endif
  long long* seed,
  double* NextStoreProfile, double* NextUpdateProfile,
  double* NextStorePos, double* NextStoreDisp
#ifdef SIGMA
  , int* sigma
#endif
  ){
  
  long i,j;
  int argctarget=0;           // Number of parameters that should be used
  char command_base[1000]=""; // string that contains the desired format of the command line
  char name[200];             // string in which the file names are written

  /*
    Format the command line and count the arguments
  */
  argctarget=0;
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]); argctarget ++;
  strcat(command_base, " ");
  strcat(command_base, "file "); argctarget ++;
  strcat(command_base, "Lx "); argctarget ++;
  strcat(command_base, "Ly "); argctarget ++;
  strcat(command_base, "N "); argctarget ++;
  strcat(command_base, "dt "); argctarget ++;
  strcat(command_base, "final_time "); argctarget ++;
#ifdef RTP
  strcat(command_base, "alpha "); argctarget ++;
#elif defined ABP
  strcat(command_base, "Dr "); argctarget ++;
#endif
#ifdef DT
  strcat(command_base, "Dt "); argctarget ++;
#endif
  strcat(command_base, "v0 "); argctarget ++;
#ifdef VOFX
  strcat(command_base, "v1 "); argctarget ++;
  strcat(command_base, "v_center1 "); argctarget ++;
  strcat(command_base, "v_center2 "); argctarget ++;
  strcat(command_base, "v_right "); argctarget ++;
  strcat(command_base, "P "); argctarget ++;
#endif
#ifndef NI
  strcat(command_base, "sigma "); argctarget ++;
  strcat(command_base, "amp "); argctarget ++;
#endif
  strcat(command_base, "seed "); argctarget ++;
  strcat(command_base, "Nbinx "); argctarget++;
  strcat(command_base, "Nbiny "); argctarget++;
  strcat(command_base, "NstepProfile "); argctarget++;
  strcat(command_base, "StoreInterProfile "); argctarget ++;
  strcat(command_base, "StoreInterPos "); argctarget ++;
  strcat(command_base, "StoreInterDisp "); argctarget ++;
#ifdef TRACK
  strcat(command_base, "Ntrack "); argctarget++;
#endif
#ifdef GIVEN_IC
  strcat(command_base, "IC_file "); argctarget++;
#endif

  strcat(command_base, "\n");

  // Check if the call to the program was correct
  if (argc != argctarget){
    printf("%s\n",command_base);
    exit(1);
  }

  // Affect variables
  i=1;
  
  // Read in file name and create files
  sprintf(name,"%s-param",argv[i]);
  output_param[0]=fopen(name,"w");
  
  sprintf(name,"%s-profile_rho",argv[i]);
  output_profile_rho[0]=fopen(name,"w");
  
  sprintf(name,"%s-profile_m",argv[i]);
  output_profile_m[0]=fopen(name,"w");

  sprintf(name,"%s-pos",argv[i]);
  output_pos[0]=fopen(name,"w");

  sprintf(name,"%s-disp",argv[i]);
  output_disp[0]=fopen(name,"w");
  
#ifdef FORCE_PROFILE
  sprintf(name,"%s-profile_force",argv[i]);
  output_profile_force[0]=fopen(name,"w");
#endif
  
#ifdef SIGMA
  sprintf(name,"%s-profile_sigma",argv[i]);
  output_profile_sigma[0]=fopen(name,"w");
#endif
  
#ifdef TRACK
  sprintf(name,"%s-traj",argv[i]);
  output_traj[0]=fopen(name,"w");
#endif
  
  printf("created files\n");
  
  i++;
  
  // Read in basic parameters
  Param[0].Lx                     = strtod(argv[i], NULL); i++;
  Param[0].Ly                     = strtod(argv[i], NULL); i++;
  Param[0].N                      = (long) strtod(argv[i], NULL); i++;
  Param[0].dt                     = strtod(argv[i], NULL); i++;
  Param[0].final_time             = strtod(argv[i], NULL); i++;
#ifdef RTP
  Param[0].alpha                  = strtod(argv[i], NULL); i++;
#elif defined ABP
  Param[0].Dr                     = strtod(argv[i], NULL); i++;
#endif
#ifdef DT
  Param[0].Dt                     = strtod(argv[i], NULL); i++;
#endif
  Param[0].v0                     = strtod(argv[i], NULL); i++;
#ifdef VOFX
  Param[0].v1                     = strtod(argv[i], NULL); i++;
  Param[0].v_center1              = strtod(argv[i], NULL); i++;
  Param[0].v_center2              = strtod(argv[i], NULL); i++;
  Param[0].v_right                = strtod(argv[i], NULL); i++;
  Param[0].P                      = (int) strtod(argv[i], NULL); i++;
#endif
#ifndef NI
  Param[0].sigma                  = strtod(argv[i], NULL); i++;
  Param[0].amp                    = strtod(argv[i], NULL); i++;
#endif
  seed[0]                         = strtod(argv[i], NULL); i++;
  Param[0].Nbinx                  = strtod(argv[i], NULL); i++;
  Param[0].Nbiny                  = strtod(argv[i], NULL); i++;
  Param[0].NstepProfile           = (int) strtod(argv[i], NULL); i++;
  Param[0].StoreInterProfile      = strtod(argv[i], NULL); i++;
  Param[0].StoreInterPos          = strtod(argv[i], NULL); i++;
  Param[0].StoreInterDisp         = strtod(argv[i], NULL); i++;
#ifdef TRACK
  Param[0].Ntrack                 = strtod(argv[i], NULL); i++;
#endif
#ifdef GIVEN_IC
  Param[0].IC_file                = argv[i]; i++;
#endif
  
#ifdef ABP
  Param[0].tau                = 1/Param[0].Dr;
#elif defined(RTP)
  Param[0].tau                = 1/Param[0].alpha;
#endif
  
  // define parameters that are functions of the inputs and/or known
  NextStoreProfile[0]      = Param[0].StoreInterProfile;
  NextUpdateProfile[0]     = NextStoreProfile[0] - Param[0].tau * (Param[0].NstepProfile-1);
  
  // Must have enough time to make enough measurements between profile storages
  assert(Param[0].NstepProfile*Param[0].tau<=Param[0].StoreInterProfile);

  NextStorePos[0]          = Param[0].StoreInterPos;
  NextStoreDisp[0]         = Param[0].StoreInterDisp;

#ifdef VOFX
  Param[0].dv              = Param[0].v1 - Param[0].v0;
  Param[0].dxR             = Param[0].v_right - Param[0].v_center2;
  Param[0].motif_L         = Param[0].Lx / ((double) Param[0].P);
  
#ifdef LANDSCAPE_CUBIC
  // use cubic interpolation from v0 to v1 and vice versa that has derivative zero at the endpoints 
  // because v_left=0, we have b1=0 and a1=v0
  Param[0].c1    = 3*Param[0].dv / pow(Param[0].v_center1,2);
  Param[0].d1    = -2*Param[0].dv / pow(Param[0].v_center1,3);
  
  Param[0].a2    = -Param[0].dv * pow(Param[0].v_center2,2) * (3*Param[0].dxR + 2*Param[0].v_center2) / pow(Param[0].dxR,3) + Param[0].v1;
  Param[0].b2    = 6*Param[0].dv * Param[0].v_center2 * (Param[0].dxR + Param[0].v_center2) / pow(Param[0].dxR,3);
  Param[0].c2    = -3*Param[0].dv * (Param[0].dxR + 2*Param[0].v_center2) / pow(Param[0].dxR,3);
  Param[0].d2    = 2*Param[0].dv / pow(Param[0].dxR,3);
#endif
#endif

#ifdef ABP
  Param[0].sqrt2Drdt       = sqrt(2*Param[0].Dr*Param[0].dt);
#endif

#ifdef DT
  Param[0].sqrt2Dtdt       = sqrt(2*Param[0].Dt*Param[0].dt);
#endif

#ifndef NI
#if defined(FORCE_HARMONIC) || defined(FORCE_INV_HARMONIC)
  Param[0].amp_over_sigma  = Param[0].amp / Param[0].sigma;
#elif defined FORCE_QUARTIC
  Param[0].force_coef1 = 4*Param[0].amp/pow(Param[0].sigma,4);
  Param[0].force_coef2 = 4*Param[0].amp/pow(Param[0].sigma,2);
#endif
#ifdef WCA
  Param[0].rmax            = Param[0].sigma * pow(2, 1/6.);
  Param[0].force_coef1     = 4*Param.amp*12*pow(Param[0].sigma,12);
  Param[0].force_coef2     = 4*Param.amp*6*pow(Param[0].sigma,6);
#else 
  Param[0].rmax            = Param[0].sigma;
#endif
  Param[0].rmax2           = Param[0].rmax*Param[0].rmax;
  Param[0].rbox            = Param[0].rmax;
  Param[0].NxBox           = (int) (floor(Param[0].Lx/Param[0].rbox)+EPS);
  Param[0].NyBox           = (int) (floor(Param[0].Ly/Param[0].rbox)+EPS);
  
  // System width must be integer multiple of interaction length
  assert((double)Param[0].NxBox - Param[0].Lx/Param[0].rbox<EPS);
  assert((double)Param[0].NyBox - Param[0].Ly/Param[0].rbox<EPS);
#endif
  
  Param[0].binwidthx       = (double) Param[0].Lx/((double) Param[0].Nbinx);
  Param[0].binwidthy       = (double) Param[0].Ly/((double) Param[0].Nbiny);
  _time[0]                 = 0;
  prev_percentage[0]       = 0;
  
  init_genrand64(seed[0]);
  
#ifdef VOFX
  // activity landscape must satisfy 0 <= v_center1 <= v_center2 <= v_right
  assert(0<=Param[0].v_center1);
  assert(Param[0].v_center1<=Param[0].v_center2);
  assert(Param[0].v_center2<=Param[0].v_right);
  
  // activity landscapes must lie inside system
  assert(Param[0].v_right*Param[0].P<=Param[0].Lx);
#endif

  // Allocate space for simulation data
  Particles[0]                = (particle*) malloc(sizeof(particle)*Param[0].N);
  Displacements[0]            = (double*) calloc(2*Param[0].N,sizeof(double));
  Integrated_Displacements[0] = (double*) calloc(2*Param[0].N,sizeof(double));
#ifndef NI
  forces[0]                   = (double*) calloc(2*Param[0].N,sizeof(double));
#endif
  
  profile_rho[0]              = (long**) malloc(Param[0].Nbinx*sizeof(long*));
  for (i=0;i<Param[0].Nbinx;i++){
    profile_rho[0][i]         = (long*) calloc(Param[0].Nbiny,sizeof(long));
    for (j=0;j<Param[0].Nbiny;j++){
      profile_rho[0][i][j] = 0;
    }
  }
  profile_m[0]              = (double**) malloc(Param[0].Nbinx*sizeof(double*));
  for (i=0;i<Param[0].Nbinx;i++){
    profile_m[0][i]         = (double*) calloc(Param[0].Nbiny,sizeof(double));
    for (j=0;j<Param[0].Nbiny;j++){
      profile_m[0][i][j] = 0.0;
    }
  }
#ifdef FORCE_PROFILE
  profile_force[0]                = (double***) malloc(Param[0].Nbinx*sizeof(double**));
  for (i=0;i<Param[0].Nbinx;i++){
    profile_force[0][i]           = (double**) calloc(Param[0].Nbiny,sizeof(double*));
    for (j=0;j<Param[0].Nbiny;j++){
      profile_force[0][i][j]      = (double*) calloc(5, sizeof(double));
      profile_force[0][i][j][0]   = 0.0;
      profile_force[0][i][j][1]   = 0.0;
      profile_force[0][i][j][2]   = 0.0;
      profile_force[0][i][j][3]   = 0.0;
      profile_force[0][i][j][4]   = 0.0;
    }
  }
#endif
  
#ifndef NI

#ifdef SIGMA
  sigma[0]                    = 1; // whether or not to update sigma at the next step
  profile_sigma[0]           = (double***) malloc(Param[0].NxBox*sizeof(double**));
  for (i=0;i<Param[0].NxBox;i++){
    profile_sigma[0][i]      = (double**) calloc(Param[0].NyBox,sizeof(double*));
    for (j=0;j<Param[0].NyBox;j++){
      profile_sigma[0][i][j] = (double*) calloc(10, sizeof(double));
      profile_sigma[0][i][j][0] = 0.0; // sigma_xx
      profile_sigma[0][i][j][1] = 0.0; // sigma_xy
      profile_sigma[0][i][j][2] = 0.0; // sigma_yx
      profile_sigma[0][i][j][3] = 0.0; // sigma_yy
      profile_sigma[0][i][j][4] = 0.0; // G_xx
      profile_sigma[0][i][j][5] = 0.0; // G_xy
      profile_sigma[0][i][j][6] = 0.0; // G_yx
      profile_sigma[0][i][j][7] = 0.0; // G_yy
      profile_sigma[0][i][j][8] = 0.0; // FA_x
      profile_sigma[0][i][j][9] = 0.0; // FA_y
    }
  }
#endif


  // Construct list of 1st particle in each box (initialized as empty)
  Boxes[0] = (long**) malloc(Param[0].NxBox*sizeof(long*));
  for (i=0;i<Param[0].NxBox;i++){
    Boxes[0][i] = (long*) calloc(Param[0].NyBox,sizeof(long));
    for (j=0;j<Param[0].NyBox;j++)
      Boxes[0][i][j]=-1;
  }
  
  // Construct empty list of neighbors
  Neighbours[0] = (long*) calloc((Param[0].N+1)*2,sizeof(long));
  for(i=0;i<Param[0].N;i++){
    Neighbours[0][2*i]=-1;
    Neighbours[0][2*i+1]=-1;
  }
  
  // Construct list of next boxes
  NeighbouringBoxes[0] = (box***) malloc(Param[0].NxBox*sizeof(box**));
  DefineNeighbouringBoxes(NeighbouringBoxes[0],Param[0]);
#endif

// initialize particle locations and polarizations  
#ifdef RANDOM_IC
  for (i=0;i<Param[0].N;i++){
  
    // random initialize positions uniformly over interval
    Particles[0][i].x=Param[0].Lx*genrand64_real2();
    Particles[0][i].y=Param[0].Ly*genrand64_real2();
    
    // randomly pick direction
    Particles[0][i].theta = 2*M_PI*genrand64_real2();
    
#ifdef RTP
    Particles[0][i].next_time = -2*log(genrand64_real3())/Param[0].alpha;
#endif
    
#ifndef NI
    // Add particle to its box and link to its neighbors
    Particles[0][i].bi= (int) (floor(Particles[0][i].x/Param[0].rbox) + EPS);
    Particles[0][i].bj= (int) (floor(Particles[0][i].y/Param[0].rbox) + EPS);
    AddinBox(i,Particles[0][i].bi,Particles[0][i].bj,Boxes[0],Neighbours[0]);
#endif
  }
  
  
#elif defined GIVEN_IC

  // open file for reading
  FILE *fp;
  fp = fopen(Param[0].IC_file, "r");

  // Read the first line to get the number of columns
  int MAX_LINE_LEN=1024;
  char line[MAX_LINE_LEN];

  // 4 columns: time, particle ID, x, theta
  int num_columns = 4;
  char* pch;

  // no header
  // Read the file, line by line
  // Note: this will throw an error if the number of lines doesn't equal N
  // It also will only assign the particle ID's that are present in column 2 (--> error later if one missing)
  double tmp, x_tmp;
  int theta_tmp;
  long id;
  char *ptr;
  while (fgets(line, MAX_LINE_LEN, fp) != NULL) {
  
    // read line to time, id, x, theta variables
    sscanf(line, "%lg\t%ld\t%lg\t%d", &tmp, &id,  &x_tmp,  &y_tmp, &theta_tmp);
    Particles[0][id].x = x_tmp;
    Particles[0][id].y = y_tmp;
    Particles[0][id].theta = theta_tmp;
    
#ifndef NI
    // Add particle to its box and link to its neighbors
    Particles[0][id].bi= (int) (floor(Particles[0][id].x/Param[0].rbox) + EPS);
    Particles[0][id].bj= (int) (floor(Particles[0][id].y/Param[0].rbox) + EPS);
    AddinBox(id,Particles[0][id].bi,Particles[0][id].bj,Boxes[0],Neighbours[0]);
#endif
  }
  
  // set the current time 
  _time[0] = tmp;
  
  NextStoreProfile[0]  += _time[0];
  NextUpdateProfile[0] += _time[0];
  NextStorePos[0]      += _time[0];
  NextStoreDisp[0]     += _time[0];

  // Close the file
  fclose(fp);
#endif
}


double vofx_function(double x, param Param){

#ifdef VOFX
  /*
    Make a piecewise linear velocity field.
    The velocity field will be repeated P times, equally spaced.
  */
  
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in
  
  n_motif = (int) (floor(x/Param.motif_L) + EPS);
  x_inside_motif = x - Param.motif_L*((double) n_motif);
  
#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  // recall that we have set v_left = 0.
  if (x_inside_motif>Param.v_right) { return Param.v0; }
  else {
    if (x_inside_motif<Param.v_center1) { return Param.v0 + Param.dv*x_inside_motif/Param.v_center1; }
    else if (x_inside_motif < Param.v_center2) { return Param.v1; }
    else { return Param.v1 - Param.dv*(x_inside_motif - Param.v_center2)/Param.dxR; }
  }
#elif defined(LANDSCAPE_CUBIC)
    // put in piecewise cubic/constant motif as if it is the entire system
  if (x_inside_motif >= Param.v_right) { return Param.v0; }
  else
    {
      if (x_inside_motif < Param.v_center1) { return Param.v0 + Param.c1*pow(x_inside_motif,2) + Param.d1*pow(x_inside_motif,3); }
      else if (x_inside_motif < Param.v_center2) { return Param.v1; }
      else { return Param.a2 + Param.b2*x_inside_motif + Param.c2*pow(x_inside_motif,2) + Param.d2*pow(x_inside_motif,3); }
    }
#endif
#else
  return Param.v0;
#endif
}


double dvofxdx_function(double x, param Param){
#ifdef VOFX
  /* 
    derivative of the piecewise linear velocity field.
  */
  
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in
  
  n_motif = (int) (floor(x/Param.motif_L) + EPS);
  x_inside_motif = x - Param.motif_L*((double) n_motif);
  
#ifdef LANDSCAPE_LINEAR
  // now, put in piecewise linear motif as if the motif is the entire system.
  // recall that we have set v_left = 0.
  if (x_inside_motif>Param.v_right) { return 0.0; }
  else {
    if (x_inside_motif<Param.v_center1) { return Param.dv/Param.v_center1; }
    else if (x_inside_motif < Param.v_center2) { return 0.0; }
    else { return -Param.dv/Param.dxR; }
  }
#elif defined(LANDSCAPE_CUBIC)
  // now, put in piecewise cubic motif as if the motif is the entire system.
  // recall that we have set v_left = 0.
  if (x_inside_motif>Param.v_right) { return 0.0; }
  else {
    if (x_inside_motif<Param.v_center1) { return 2*Param.c1*x_inside_motif + 3*Param.d1*pow(x_inside_motif,2); }
    else if (x_inside_motif < Param.v_center2) { return 0.0; }
    else { return Param.b2 + 2*Param.c2*x_inside_motif + 3*Param.d2*pow(x_inside_motif,2); }
  }
#endif
#else
  return 0.0;
#endif

}




#ifndef NI

// Force_harmonic computes the repulsive harmonic force exerted by the particle k onto particle j
// with a concave-down harmonic potential
double Force_harmonic(double dist, param Param){
#ifdef FORCE_HARMONIC
  return Param.amp_over_sigma*dist;
#endif
}

// Force_inv_harmonic computes the repulsive harmonic force exerted by the particle k onto particle j
// with a concave-up harmonic potential
double Force_inv_harmonic(double dist, param Param){
#ifdef FORCE_INV_HARMONIC
  return Param.amp - Param.amp_over_sigma * dist;
#endif
}

double dfxdx_inv_harmonic(double dx2, double dist, double dist3, param Param){
#ifdef FORCE_INV_HARMONIC
  return -Param.amp_over_sigma + Param.amp/dist - Param.amp*dx2/dist3;
#endif
}

double dfxdy_inv_harmonic(double dx, double dy, double dist3, param Param){
#ifdef FORCE_INV_HARMONIC
  return -Param.amp*dx*dy/dist3;
#endif
}

double Force_quartic(double dx, double dx2, param Param){
#ifdef FORCE_QUARTIC
  return -Param.force_coef1*dx2*dx + Param.force_coef2*dx;
#endif
}

double dfxdx_quartic(double dx, double dx2, double dist, double dist2, param Param){
#ifdef FORCE_QUARTIC
  return -Param.force_coef1*dist2 - Param.force_coef1*2*dx2 + Param.force_coef2;
#endif
}
double dfxdy_quartic(double dx, double dy, param Param){
#ifdef FORCE_QUARTIC
  return -2*Param.force_coef1*dx*dy;
#endif
}


//compute the repulsive force due to a piecewise linear "tent" interaction potential
double Force_linear(param Param){
#ifdef FORCE_LINEAR
  return Param.amp;
#endif
}

//compute the force due to WCA potential (cutoff LJ)
//force_coef1 = 4*amp*12*sigma^12 
//force_coef2 = 4*amp*6*sigma^6
void Force_WCA(double dist, double dist2,param Param){
#ifdef WCA
  double dist7=dist2*dist2*dist2*dist;
  double dist13=dist7*dist2*dist2*dist2;
  return Param.force_coef1/dist13 - Param.force_coef2/dist7;
#endif
}


//This function computes the interactions between particle J and all
//the subsequent particles in the same box
void Compute_force_same_box(long j,long* Neighbours,double* forces, param Param, particle* Particles
#ifdef SIGMA
                              ,int sigma, double**** profile_sigma
#endif
){
  long k; // index of neighbors
  double Force_x, Force_y; // Force exerted between two particles
  double dx,dy,dx2,dy2; // displacement between particles
  double dist; // distance between particles 
  double dist2,dist3; // distance between particles squared
  double force;   // distance between particles squared
#ifdef SIGMA
  int bi, bj;     // box index
  double sigma_xx, sigma_xy, sigma_yx, sigma_yy; // stress tensor components
#endif
  
  // Loop through all the following neighbors of j
  k=Neighbours[2*j+1]; // start with the 1st neighbor
  
  // As long as there are neighbors, iterate
  while(k!=-1){
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
    dx2 = dx*dx;
    dy2 = dy*dy;
    dist2 = dx2 + dy2;
    
    //If the particles are closer than Param.rmax, they interact
    if (dist2<Param.rmax2){
#ifdef FORCE_HARMONIC
      dist = pow(dist2,0.5);
      force=Force_harmonic(dist, Param);
#elif defined FORCE_INV_HARMONIC
      dist = pow(dist2,0.5);
      force=Force_inv_harmonic(dist, Param);
#elif defined FORCE_QUARTIC
      dist = pow(dist2,0.5);
      force=Force_quartic(dist, dist2, Param);
#elif defined FORCE_LINEAR
      force=Force_linear(Param);
#elif defined WCA
      dist = pow(dist2,0.5);
      force=Force_WCA(dist, dist2, Param);
#endif
      Force_x = force*dx/dist;
      Force_y = force*dy/dist;
      forces[2*j] += Force_x;
      forces[2*j+1] += Force_y;
      forces[2*k] -= Force_x;
      forces[2*k+1] -= Force_y;
    
#ifdef SIGMA
      // measure stress tensor
      if (sigma){
        sigma_xx = dx*Force_x;
        sigma_xy = dx*Force_y;
        sigma_yx = sigma_xy;
        sigma_yy = dy*Force_y;
      
        bi = (int) (floor(Particles[j].x/Param.rbox)+EPS);
        bj = (int) (floor(Particles[j].y/Param.rbox)+EPS);
      
        profile_sigma[0][bi][bj][0] += sigma_xx;
        profile_sigma[0][bi][bj][1] += sigma_xy;
        profile_sigma[0][bi][bj][2] += sigma_yx;
        profile_sigma[0][bi][bj][3] += sigma_yy;
        
        // add active part in update_particles function
      }
#endif
    }
    k=Neighbours[2*k+1]; // k is now the next particle in the list
  }
}

// This functions compute the interactions between particle j and all
// the particles in the box situated on its right, using periodic
// boundary conditions
void Compute_force_neighbours(long j,int bi,int bj,box*** NeighbouringBoxes,particle* Particles,
                              param Param,double* forces, long* Neighbours, long** Boxes
#ifdef SIGMA
                              ,int sigma, double**** profile_sigma
#endif
#ifdef EPR
                              , int mid,double**** dists_mid, double**** fps_mid, int*** interact_mid
#endif
){
  int nbi,nbj;     //box name
  int m;           // neighbour box id
  double dx_box;   // particle offset
  double dy_box;   // particle offset
  long k;          //particle index
  double Force_x,Force_y,force;    //Force between particles
  double dx,dy,dx2,dy2,dist,dist2,dist3;       // distance between particles
  
#ifdef SIGMA
  double sigma_xx, sigma_xy, sigma_yx, sigma_yy; // stress tensor components
  double frac_0, frac_2, frac_4; // frac_k = fraction of separation vector in box k
  double sigma_xx_part, sigma_xy_part, sigma_yx_part, sigma_yy_part; // stress tensor components (piece within box)
  double sigma_xx_part1, sigma_xy_part1, sigma_yx_part1, sigma_yy_part1; // stress tensor components (piece within neighbor box)
  double intersect_x; // where it intersects the edge of the box (x coordinate)
  double box_left, box_bottom, box_right, box_top; // coordinates of the box sides
  int nbi1,nbj1; // overlapping box
#endif
  
  // NeighbouringBoxes[bi][bj][0] is box (bi,bj) itself
  for (m=1;m<5;m++){
    nbi=NeighbouringBoxes[bi][bj][m].i;        // x index of the mth neighbour
    nbj=NeighbouringBoxes[bi][bj][m].j;        // y index of the mth neighbour
    dx_box=NeighbouringBoxes[bi][bj][m].epsilonx;// x offset to be added to the particles in
                                              // box (bi,bj) to take into account the periodic
                                              // boundary conditions
    dy_box=NeighbouringBoxes[bi][bj][m].epsilony;// y offset to be added to the particles in
                                  // box (bi,bj) to take into account the periodic
                                  // boundary conditions
                                  
#ifdef SIGMA
    // get coordinates of edges of box for geometric calculations below
    if (sigma){
      box_left = Param.rbox * (floor(Particles[j].x/Param.rbox)+EPS);
      box_bottom = Param.rbox * (floor(Particles[j].y/Param.rbox)+EPS);
      box_right = box_left + Param.rbox;
      box_top = box_bottom + Param.rbox;
    }
#endif

    //Loop through all particles k in box n
    k=Boxes[nbi][nbj]; //Start with the first
  
    //As long as there are particles in the box, iterate
    while(k!=-1){

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
      dx2=dx*dx;
      dy2=dy*dy;
      dist2=dx2+dy2;
    
      //If the particles are closer than Param.rmax, they interact
      if (dist2<Param.rmax2){
#ifdef FORCE_HARMONIC
        dist=pow(dist2,0.5);
        force=Force_harmonic(dist, Param);
#elif defined FORCE_INV_HARMONIC
        dist=pow(dist2,0.5);
        force=Force_inv_harmonic(dist, Param);
#elif defined FORCE_QUARTIC
        dist = pow(dist2,0.5);
        force=Force_quartic(dist, dist2, Param);
#elif defined FORCE_LINEAR
        force=Force_linear(Param);
#elif defined WCA
        dist=pow(dist2,0.5);
        force=Force_WCA(dist, dist2, Param);
#endif
        Force_x = force*dx/dist;
        Force_y = force*dy/dist;
        forces[2*j] += Force_x;
        forces[2*j+1] += Force_y;
        forces[2*k] -= Force_x;
        forces[2*k+1] -= Force_y;
        
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
            sigma_yx_part = sigma_xy_part;
            sigma_yy_part = sigma_yy*frac_0;
            profile_sigma[0][bi][bj][0] += sigma_xx_part;
            profile_sigma[0][bi][bj][1] += sigma_xy_part;
            profile_sigma[0][bi][bj][2] += sigma_yx_part;
            profile_sigma[0][bi][bj][3] += sigma_yy_part;
            
            // add to neighbor box
            profile_sigma[0][nbi][nbj][0] += sigma_xx-sigma_xx_part;
            profile_sigma[0][nbi][nbj][1] += sigma_xy-sigma_xy_part;
            profile_sigma[0][nbi][nbj][2] += sigma_yx-sigma_yx_part;
            profile_sigma[0][nbi][nbj][3] += sigma_yy-sigma_yy_part;
            
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
            sigma_yx_part = sigma_xy_part;
            sigma_yy_part = sigma_yy*frac_0;
            profile_sigma[0][bi][bj][0] += sigma_xx_part;
            profile_sigma[0][bi][bj][1] += sigma_xy_part;
            profile_sigma[0][bi][bj][2] += sigma_yx_part;
            profile_sigma[0][bi][bj][3] += sigma_yy_part;
            
            // add to neighbor box
            sigma_xx_part1 = sigma_xx*frac_2;
            sigma_xy_part1 = sigma_xy*frac_2;
            sigma_yx_part1 = sigma_xy_part1;
            sigma_yy_part1 = sigma_yy*frac_2;
            profile_sigma[0][nbi][nbj][0] += sigma_xx_part1;
            profile_sigma[0][nbi][nbj][1] += sigma_xy_part1;
            profile_sigma[0][nbi][nbj][2] += sigma_yx_part1;
            profile_sigma[0][nbi][nbj][3] += sigma_yy_part1;
            
            // add to overlapping box
            profile_sigma[0][nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
            profile_sigma[0][nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
            profile_sigma[0][nbi1][nbj1][2] += sigma_yx-sigma_yx_part-sigma_yx_part1;
            profile_sigma[0][nbi1][nbj1][3] += sigma_yy-sigma_yy_part-sigma_yy_part1;
            
          } else if (m==3){ //to the right
            // dx will be negative
            frac_0 = (box_right-Particles[j].x)/(-dx);
            
            // add to same box
            sigma_xx_part = sigma_xx*frac_0;
            sigma_xy_part = sigma_xy*frac_0;
            sigma_yx_part = sigma_xy_part;
            sigma_yy_part = sigma_yy*frac_0;
            profile_sigma[0][bi][bj][0] += sigma_xx_part;
            profile_sigma[0][bi][bj][1] += sigma_xy_part;
            profile_sigma[0][bi][bj][2] += sigma_yx_part;
            profile_sigma[0][bi][bj][3] += sigma_yy_part;
            
            // add to neighbor box
            profile_sigma[0][nbi][nbj][0] += sigma_xx-sigma_xx_part;
            profile_sigma[0][nbi][nbj][1] += sigma_xy-sigma_xy_part;
            profile_sigma[0][nbi][nbj][2] += sigma_yx-sigma_yx_part;
            profile_sigma[0][nbi][nbj][3] += sigma_yy-sigma_yy_part;
            
          } else if (m==4){ //below and to the right
            // dx will be negative and dy will be positive
            intersect_x = (Particles[j].y - box_bottom)*(-dx)/dy;
            if (intersect_x < box_right) { // goes thru box "5"
              frac_0 = (Particles[j].y - box_bottom)/dy;
              frac_4 = (Particles[k].x+dx_box - box_right)/(-dx);
              
              // will overlap the box below.
              nbi1 = bi;
              nbj1 = (bj-1+Param.NyBox)%Param.NyBox;
              
            } else { // goes thru box 3
              frac_0 = (box_right - Particles[j].x)/(-dx);
              frac_4 = (box_bottom - (Particles[k].y+dy_box))/dy;
              
              nbi1 = NeighbouringBoxes[bi][bj][3].i;
              nbj1 = NeighbouringBoxes[bi][bj][3].j;
            }
            
            // add to same box
            sigma_xx_part = sigma_xx*frac_0;
            sigma_xy_part = sigma_xy*frac_0;
            sigma_yx_part = sigma_xy_part;
            sigma_yy_part = sigma_yy*frac_0;
            profile_sigma[0][bi][bj][0] += sigma_xx_part;
            profile_sigma[0][bi][bj][1] += sigma_xy_part;
            profile_sigma[0][bi][bj][2] += sigma_yx_part;
            profile_sigma[0][bi][bj][3] += sigma_yy_part;
            
            // add to neighbor box
            sigma_xx_part1 = sigma_xx*frac_4;
            sigma_xy_part1 = sigma_xy*frac_4;
            sigma_yx_part1 = sigma_xy_part1;
            sigma_yy_part1 = sigma_yy*frac_4;
            profile_sigma[0][nbi][nbj][0] += sigma_xx_part1;
            profile_sigma[0][nbi][nbj][1] += sigma_xy_part1;
            profile_sigma[0][nbi][nbj][2] += sigma_yx_part1;
            profile_sigma[0][nbi][nbj][3] += sigma_yy_part1;
            
            // add to overlapping box
            profile_sigma[0][nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
            profile_sigma[0][nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
            profile_sigma[0][nbi1][nbj1][2] += sigma_yx-sigma_yx_part-sigma_yx_part1;
            profile_sigma[0][nbi1][nbj1][3] += sigma_yy-sigma_yy_part-sigma_yy_part1;
            
          }
        }
#endif
      }
      k=Neighbours[2*k+1]; // k is now the next particle in the list
    }
  }
}

void Loop_Force(particle* Particles, param Param, double* forces,long* Neighbours, 
                box*** NeighbouringBoxes, long** Boxes
#ifdef SIGMA
                ,int sigma, double**** profile_sigma
#endif
#ifdef EPR
                 , int mid,double**** dists_mid, double**** fps_mid, int*** interact_mid
#endif
){
  int bi,bj; // Box index
  long j; // Particle index

  //Initialize the array of forces to zero
  memset(forces,0,2*Param.N*sizeof(double));
  
  /* 
     To compute the force on each particle, we loop through all the boxes and compute:
     - the interactions between particles inside the box
     - the interactions between particles inside the box and inside a neighboring box
  */

  // Loop through all boxes
  for(bi=0; bi<Param.NxBox; bi++){
    for(bj=0; bj<Param.NyBox; bj++){
  
      //Compute the force inside the box
      //Loop through all the particles 'j' in the box.
      j=Boxes[bi][bj]; //Start with j being the 1st particle

      //As long as j is not -1, compute its interactions with all the other particles
      while(j!=-1){
      
        // Loop first with the particles in the same box
        Compute_force_same_box(j,Neighbours,forces,Param,Particles
#ifdef SIGMA
                              ,sigma, profile_sigma
#endif
                              );

        // Loop through the neigboring box
        Compute_force_neighbours(j, bi, bj, NeighbouringBoxes, Particles, Param, forces, Neighbours, Boxes
#ifdef SIGMA
                                 ,sigma, profile_sigma
#endif
                                );
        
        //We have included all the interactions between j and other particles. Now iterate over j
        j=Neighbours[2*j+1];
      }
    // We are done with all particles in Box[i], iterate over the box
    }
  }
}


void Update_Particles(param Param,double* Displacements,particle* Particles, 
                      double* Integrated_Displacements, double* forces, long* Neighbours, 
                      box*** NeighbouringBoxes, long** Boxes,
#ifdef SIGMA
                      int sigma, double**** profile_sigma,
#endif
                      double _time){
  long n;
  int newbi,newbj;
  double v, cos_th, sin_th;
#ifdef SIGMA
  int bi,bj;
  double vp, v_Dr_dt, vp_Dr_dt;
#endif
  
  // incorporate interactions with spatial hashing
  Loop_Force(Particles, Param, forces, Neighbours, NeighbouringBoxes, Boxes
#ifdef SIGMA
                              ,sigma, profile_sigma
#endif
                                );
  
  for ( n=0 ; n<Param.N ; n++ ){
  
#ifdef RTP
    // If the particle tumbles, it acquires a new direction
    if ( Particles[n].next_time < _time + Param.dt){
      Particles[n].theta = 2*M_PI*genrand64_real3();
      Particles[n].next_time += -2*log(genrand64_real3())/Param.alpha;
    }
#elif defined ABP
    Particles[n].theta += Param.sqrt2Drdt * gasdev();
    
    if(Particles[n].theta>2*M_PI) Particles[n].theta-=2*M_PI;
    if(Particles[n].theta<0) Particles[n].theta+=2*M_PI;
#endif

    // Compute the displacement of each particle as x += v cos(theta) dt, y += v sin(theta) dt
    v = vofx_function(Particles[n].x, Param);
    cos_th = cos(Particles[n].theta);
    sin_th = sin(Particles[n].theta);
    Displacements[2*n]   = v * cos_th * Param.dt;
    Displacements[2*n+1] = v * sin_th * Param.dt;

#ifdef DT
    Displacements[2*n]   += Param.sqrt2Dtdt * gasdev();
    Displacements[2*n+1] += Param.sqrt2Dtdt * gasdev();
#endif

    // incorporate interactions
    Displacements[2*n]   += forces[2*n]*Param.dt;
    Displacements[2*n+1] += forces[2*n+1]*Param.dt;
    
#ifdef SIGMA
    if (sigma){
      v = vofx_function(Particles[n].x, Param);
      vp = dvofxdx_function(Particles[n].x, Param);
      bi = Particles[n].bi;
      bj = Particles[n].bj;
      v_Dr_dt = (v*Param.tau)/Param.dt;
      vp_Dr_dt = (vp*Param.tau)/Param.dt;
      
      // add active stress 
      profile_sigma[0][bi][bj][4] += cos_th * Displacements[2*n] * v_Dr_dt; // G_xx
      profile_sigma[0][bi][bj][5] += cos_th * Displacements[2*n+1] * v_Dr_dt; // G_xy
      profile_sigma[0][bi][bj][6] += sin_th * Displacements[2*n] * v_Dr_dt; // G_yx
      profile_sigma[0][bi][bj][7] += sin_th * Displacements[2*n+1] * v_Dr_dt; // G_yy
      
      // add active force 
      profile_sigma[0][bi][bj][8] += cos_th * Displacements[2*n] * vp_Dr_dt; // FA_x
      profile_sigma[0][bi][bj][9] += sin_th * Displacements[2*n] * vp_Dr_dt; // FA_y
    }
#endif
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[2*n];
    Particles[n].y += Displacements[2*n+1];

    // Take care of periodic boundary conditions
#ifdef PBC
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;
#endif

    // Update box membership
    newbi = (int) (floor(Particles[n].x/Param.rbox) + EPS);
    newbj = (int) (floor(Particles[n].y/Param.rbox) + EPS);

    if(Particles[n].bi!=newbi || Particles[n].bj!=newbj){
      RemovefromBox(n,Particles[n].bi,Particles[n].bj,Boxes,Neighbours);
      AddinBox(n,newbi,newbj,Boxes,Neighbours);
      Particles[n].bi=newbi;
      Particles[n].bj=newbj;
    }

    // Update the cumulated displacements
    Integrated_Displacements[2*n] += Displacements[2*n];
    Integrated_Displacements[2*n+1] += Displacements[2*n+1];
  }
}
#endif


#ifdef NI
void Update_Particles_NI(param Param,double* Displacements,particle* Particles, double* Integrated_Displacements, double _time){
  long n;
  double v;
  
  for ( n=0 ; n<Param.N ; n++ ){
#ifdef RTP
    // If the particle tumbles, it acquires a new direction
    if ( Particles[n].next_time < _time + Param.dt){
      Particles[n].theta = 2*M_PI*genrand64_real3();
      Particles[n].next_time += -2*log(genrand64_real3())/Param.alpha;
    }
#elif defined ABP
    Particles[n].theta += Param.sqrt2Drdt * gasdev();
    
    if(Particles[n].theta>2*M_PI) Particles[n].theta-=2*M_PI;
    if(Particles[n].theta<0) Particles[n].theta+=2*M_PI;
#endif

    // Compute the displacement of each particle as x += v theta dt where theta= +/- 1
    v = vofx_function(Particles[n].x, Param);
    Displacements[2*n]   = v * cos(Particles[n].theta) * Param.dt;
    Displacements[2*n+1] = v * sin(Particles[n].theta) * Param.dt;

#ifdef DT
    Displacements[2*n]   += Param.sqrt2Dtdt * gasdev();
    Displacements[2*n+1] += Param.sqrt2Dtdt * gasdev();
#endif

  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[2*n];
    Particles[n].y += Displacements[2*n+1];

    // Take care of periodic boundary conditions
#ifdef PBC
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;
#endif

    // Update the cumulated displacements
    Integrated_Displacements[2*n] += Displacements[2*n];
    Integrated_Displacements[2*n+1] += Displacements[2*n+1];
  }
}
#endif


/*
	Print percentage of completed simulation
*/
void PrintSimulationProgress(double _time,double FinalTime,double * prev_percentage){
  int progress;
  
  progress = (int) (_time*100/FinalTime);
  if (progress>=prev_percentage[0]+5) {
    printf("\rIn progress [%d %%]", progress);
    fflush(stdout);
    prev_percentage[0]=progress;
  }
}

void Store_Parameters(int argc, char* argv[], FILE* output_param, param Param,
  long long seed){
  long i;
  
  /* Store parameters */
  for(i=0;i<argc;i++){
    fprintf(output_param,"%s ",argv[i]);
  }
  
  fprintf(output_param,"\n*In this version, profile update intervals was changed to single persistence time\n");
  
  fprintf(output_param,"Lx is %lg\n", Param.Lx);
  fprintf(output_param,"Ly is %lg\n", Param.Ly);
  fprintf(output_param,"N is %ld\n", Param.N);
  fprintf(output_param,"dt is %lg\n", Param.dt);
  fprintf(output_param,"final_time is %lg\n", Param.final_time);
#ifdef RTP
  fprintf(output_param,"alpha is %lg\n", Param.alpha);
#elif defined ABP
  fprintf(output_param,"Dr is %lg\n", Param.Dr);
#endif
#ifdef DT
  fprintf(output_param,"Dt is %lg\n", Param.Dt);
#endif
  fprintf(output_param,"v0 is %lg\n", Param.v0);
#ifdef VOFX
  fprintf(output_param,"v1 is %lg\n", Param.v1);
  fprintf(output_param,"v_center1 is %lg\n", Param.v_center1);
  fprintf(output_param,"v_center2 is %lg\n", Param.v_center2);
  fprintf(output_param,"v_right is %lg\n", Param.v_right);
  fprintf(output_param,"P is %d\n", Param.P);
#endif
#ifndef NI
  fprintf(output_param,"sigma is %lg\n", Param.sigma);
  fprintf(output_param,"amp is %lg\n", Param.amp);
#endif
  fprintf(output_param,"seed is %lld\n", seed);
  fprintf(output_param,"Nbinx is %d\n", Param.Nbinx);
  fprintf(output_param,"Nbiny is %d\n", Param.Nbiny);
  fprintf(output_param,"NstepProfile is %d\n", Param.NstepProfile);
  fprintf(output_param,"StoreInterProfile is %lg\n", Param.StoreInterProfile);
  fprintf(output_param,"StoreInterPos is %lg\n", Param.StoreInterPos);
  fprintf(output_param,"StoreInterDisp is %lg\n", Param.StoreInterDisp);
#ifdef TRACK
  fprintf(output_param,"Ntrack is %lg\n", Param.Ntrack);
#endif
#ifdef GIVEN_IC
  fprintf(output_param,"IC_file is %lg\n", Param.IC_file);
#endif
  fflush(output_param);
}


// Store the identities, positions, and orientations of all particles, with the time
void Store_Pos(param Param,particle* Particles,double _time,FILE* output_pos){

  long n;
  for (n=0;n<Param.N;n++){
    fprintf(output_pos,"%lg\t%ld\t%lg\t%lg\t%lg\n",_time,n,Particles[n].x,Particles[n].y,Particles[n].theta);
  }
  fflush(output_pos);
}

// Store the identities and integrated displacements of all particles, with the time
void Store_Disp(param Param, double* Integrated_Displacements, double _time,FILE* output_disp){

  long n;
  for (n=0;n<Param.N;n++){
    fprintf(output_disp,"%lg\t%ld\t%lg\t%lg\n",_time,n,Integrated_Displacements[2*n],Integrated_Displacements[2*n+1]);
  }
  fflush(output_disp);
}

// Update density and magnetization profiles (add one to bin for each particle currently in it)
void Update_Profile(param Param, particle* Particles, long*** profile_rho, double*** profile_m){
  long i; // particle index
  int j,k;  // index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    profile_rho[0][j][k] += 1;
    profile_m[0][j][k] += cos(Particles[i].theta); // save x magnetization
  }
}

// Update profile of active force , other observables that add to make J
#ifdef FORCE_PROFILE
void Update_Force_Profile(param Param, particle* Particles, double* Displacements, double**** profile_force
#ifndef NI
    , double* forces
#endif
    ){
  long i;
  int j,k;
  double v, vp, fx, cos_theta;
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    
    v = vofx_function(Particles[i].x, Param);
    vp= dvofxdx_function(Particles[i].x, Param);
    cos_theta = cos(Particles[i].theta);
    
    // active force 
    profile_force[0][j][k][0] += v*cos_theta;
    
    // v(x_i) v'(x_i) cos^2 (theta_i) / D_r
    profile_force[0][j][k][1] += v*vp*pow(cos_theta, 2)*Param.tau;
    
    // (v/D_r) \dot{r}_i^\beta cos(theta)
    profile_force[0][j][k][3] += v*cos_theta*(Displacements[2*i]/Param.dt)*Param.tau;
    profile_force[0][j][k][4] += v*cos_theta*(Displacements[2*i+1]/Param.dt)*Param.tau;
    
    // v'(x_i) cos(theta_i) \sum_{j\neq i} f^x(r_i-r_j) / Dr
#ifndef NI
    fx = forces[2*i]; // x force on particle due to all other particles
    profile_force[0][j][k][2] += fx*vp*cos_theta*Param.tau;
    
    // gradient of Irving-Kirkwood stress tensor
    profile_force[0][j][k][5] += fx;
#endif
  }
}
#endif


#ifdef SIGMA
void Store_Sigma(param Param, double**** profile_sigma, double _time, FILE* output_profile_sigma) {
  int j,k;
  
  for (j=0; j<Param.NxBox; j++) {
    for (k=0; k<Param.NyBox; k++) {
      fprintf(output_profile_sigma, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", _time, ((double)j)*Param.rbox, ((double)k)*Param.rbox, profile_sigma[0][j][k][0], profile_sigma[0][j][k][1], profile_sigma[0][j][k][2], profile_sigma[0][j][k][3], profile_sigma[0][j][k][4], profile_sigma[0][j][k][5], profile_sigma[0][j][k][6], profile_sigma[0][j][k][7], profile_sigma[0][j][k][8], profile_sigma[0][j][k][9]);
      profile_sigma[0][j][k][0]=0.0;
      profile_sigma[0][j][k][1]=0.0;
      profile_sigma[0][j][k][2]=0.0;
      profile_sigma[0][j][k][3]=0.0;
      profile_sigma[0][j][k][4]=0.0;
      profile_sigma[0][j][k][5]=0.0;
      profile_sigma[0][j][k][6]=0.0;
      profile_sigma[0][j][k][7]=0.0;
      profile_sigma[0][j][k][8]=0.0;
      profile_sigma[0][j][k][9]=0.0;
    }
  }
  fflush(output_profile_sigma);
}
#endif


// Store the profile of positions and polarizations
void Store_Profile(param Param, long*** profile_rho, double*** profile_m, double _time, FILE* output_profile_rho, FILE* output_profile_m){
  int j,k; // profile index
  
  // columns: time, x, y, count
  for (j=0; j<Param.Nbinx; j++){
    for (k=0; k<Param.Nbiny; k++){
      fprintf(output_profile_rho, "%lg\t%lg\t%lg\t%ld\n", _time, ((double)j)*Param.binwidthx, ((double)k)*Param.binwidthy, profile_rho[0][j][k]);
      profile_rho[0][j][k]=0;
      fprintf(output_profile_m, "%lg\t%lg\t%lg\t%lg\n", _time, ((double)j)*Param.binwidthx, ((double)k)*Param.binwidthy, profile_m[0][j][k]);
      profile_m[0][j][k]=0.0;
    }
  }
  fflush(output_profile_rho);
  fflush(output_profile_m);
}


#ifdef FORCE_PROFILE
void Store_Force_Profile(param Param, double**** profile_force, double _time, FILE* output_profile_force){
  int j,k; // profile index
  
  // columns: time, x, y, count
  for (j=0; j<Param.Nbinx; j++){
    for (k=0; k<Param.Nbiny; k++){
      fprintf(output_profile_force, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", 
              _time, ((double)j)*Param.binwidthx, ((double)k)*Param.binwidthy, 
              profile_force[0][j][k][0], profile_force[0][j][k][1], 
              profile_force[0][j][k][2], profile_force[0][j][k][3], 
              profile_force[0][j][k][4], profile_force[0][j][k][5]);
      profile_force[0][j][k][0]=0.0;
      profile_force[0][j][k][1]=0.0;
      profile_force[0][j][k][2]=0.0;
      profile_force[0][j][k][3]=0.0;
      profile_force[0][j][k][4]=0.0;
      profile_force[0][j][k][5]=0.0;
    }
  }
  fflush(output_profile_force);
}
#endif

#ifdef TRACK
// Store the complete trajectories (position and polarization) of first Ntrack particles
void Store_Trajectories(param Param, particle* Particles, FILE* output_traj, double _time){
  long n; //particle index
  
  for (n=0; n<Param.Ntrack; n++){
    fprintf(output_traj,"%lg\t%ld\t%lg\t%lg\t%lg\n",_time,n,Particles[n].x,Particles[n].y,Particles[n].theta);
  }
}
#endif

