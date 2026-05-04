

#ifndef NI
#include "spatial_hashing.c"
#endif


void Initialize_parameters(int argc, char* argv[],FILE** output_param, FILE** output_profile_rho, 
  FILE** output_pos, FILE** output_disp, param* Param, double* _time, double* prev_percentage,
  particle** Particles,double** Displacements, double** Integrated_Displacements, 
#ifndef NI
  double** forces, long*** Boxes, long** Neighbours, box**** NeighbouringBoxes,
#endif
  long*** profile_rho,long long* seed,double* NextStoreProfile, double* NextUpdateProfile,
  double* NextStorePos, double* NextStoreDisp
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
  strcat(command_base, "v0 "); argctarget ++;
  strcat(command_base, "v1 "); argctarget ++;
  strcat(command_base, "v_center1 "); argctarget ++;
  strcat(command_base, "v_center2 "); argctarget ++;
  strcat(command_base, "v_right "); argctarget ++;
  strcat(command_base, "P "); argctarget ++;
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

  sprintf(name,"%s-pos",argv[i]);
  output_pos[0]=fopen(name,"w");

  sprintf(name,"%s-disp",argv[i]);
  output_disp[0]=fopen(name,"w");
  
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
  Param[0].v0                     = strtod(argv[i], NULL); i++;
  Param[0].v1                     = strtod(argv[i], NULL); i++;
  Param[0].v_center1              = strtod(argv[i], NULL); i++;
  Param[0].v_center2              = strtod(argv[i], NULL); i++;
  Param[0].v_right                = strtod(argv[i], NULL); i++;
  Param[0].P                      = (int) strtod(argv[i], NULL); i++;
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

  Param[0].dv              = Param[0].v1 - Param[0].v0;
  Param[0].dxR             = Param[0].v_right - Param[0].v_center2;
  Param[0].motif_L         = Param[0].Lx / ((double) Param[0].P);
  
  // use cubic interpolation from v0 to v1 and vice versa that has derivative zero at the endpoints 
  // because v_left=0, we have b1=0 and a1=v0
  Param[0].c1    = 3*Param[0].dv / pow(Param[0].v_center1,2);
  Param[0].d1    = -2*Param[0].dv / pow(Param[0].v_center1,3);
  
  Param[0].a2    = -Param[0].dv * pow(Param[0].v_center2,2) * (3*Param[0].dxR + 2*Param[0].v_center2) / pow(Param[0].dxR,3) + Param[0].v1;
  Param[0].b2    = 6*Param[0].dv * Param[0].v_center2 * (Param[0].dxR + Param[0].v_center2) / pow(Param[0].dxR,3);
  Param[0].c2    = -3*Param[0].dv * (Param[0].dxR + 2*Param[0].v_center2) / pow(Param[0].dxR,3);
  Param[0].d2    = 2*Param[0].dv / pow(Param[0].dxR,3);


#ifdef ABP
  Param[0].sqrt2Drdt       = sqrt(2*Param[0].Dr*Param[0].dt);
#endif

#ifndef NI
  Param[0].amp_over_sigma  = Param[0].amp / Param[0].sigma;
  Param[0].rmax            = Param[0].sigma;
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
  
  // activity landscape must satisfy 0 <= v_center1 <= v_center2 <= v_right
  assert(0<=Param[0].v_center1);
  assert(Param[0].v_center1<=Param[0].v_center2);
  assert(Param[0].v_center2<=Param[0].v_right);
  
  // activity landscapes must lie inside system
  assert(Param[0].v_right*Param[0].P<=Param[0].Lx);

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
  
#ifndef NI

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
  DefineNeighbouringBoxes(NeighbouringBoxes[0],Param[0].NxBox, Param[0].NyBox, Param[0].Lx, Param[0].Ly);
#endif

// initialize particle locations and polarizations  
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
}


double vofx_function(double x, param Param){
  /*
    Make a piecewise linear velocity field.
    The velocity field will be repeated P times, equally spaced.
  */
  
  int n_motif;           // which motif it lies in
  double x_inside_motif; // position within the motif it lies in
  
  n_motif = (int) (floor(x/Param.motif_L) + EPS);
  x_inside_motif = x - Param.motif_L*((double) n_motif);
  
    // put in piecewise cubic/constant motif as if it is the entire system
  if (x_inside_motif >= Param.v_right) { return Param.v0; }
  else
    {
      if (x_inside_motif < Param.v_center1) { return Param.v0 + Param.c1*pow(x_inside_motif,2) + Param.d1*pow(x_inside_motif,3); }
      else if (x_inside_motif < Param.v_center2) { return Param.v1; }
      else { return Param.a2 + Param.b2*x_inside_motif + Param.c2*pow(x_inside_motif,2) + Param.d2*pow(x_inside_motif,3); }
    }
}




#ifndef NI

// Force_harmonic computes the repulsive harmonic force exerted by the particle k onto particle j
// with a concave-down harmonic potential
double Force_harmonic(double dist, param Param){
  return Param.amp_over_sigma*dist;
}


//This function computes the interactions between particle J and all
//the subsequent particles in the same box
void Compute_force_same_box(long j,long* Neighbours,double* forces, param Param, particle* Particles){
  long k; // index of neighbors
  double Force_x, Force_y; // Force exerted between two particles
  double dx,dy,dx2,dy2; // displacement between particles
  double dist; // distance between particles 
  double dist2,dist3; // distance between particles squared
  double force;   // distance between particles squared
  
  // Loop through all the following neighbors of j
  k=Neighbours[2*j+1]; // start with the 1st neighbor
  
  // As long as there are neighbors, iterate
  while(k!=-1){
    dx = Particles[j].x - Particles[k].x; // positive if x_j>x_k
    dy = Particles[j].y - Particles[k].y; // positive if y_j>y_k
    dx2 = dx*dx;
    dy2 = dy*dy;
    dist2 = dx2 + dy2;
    
    //If the particles are closer than Param.rmax, they interact
    if (dist2<Param.rmax2){
      dist = pow(dist2,0.5);
      force=Force_harmonic(dist, Param);
      Force_x = force*dx/dist;
      Force_y = force*dy/dist;
      forces[2*j] += Force_x;
      forces[2*j+1] += Force_y;
      forces[2*k] -= Force_x;
      forces[2*k+1] -= Force_y;
    }
    k=Neighbours[2*k+1]; // k is now the next particle in the list
  }
}

// This functions compute the interactions between particle j and all
// the particles in the box situated on its right, using periodic
// boundary conditions
void Compute_force_neighbours(long j,int bi,int bj,box*** NeighbouringBoxes,particle* Particles,
                              param Param,double* forces, long* Neighbours, long** Boxes){
  int nbi,nbj;     //box name
  int m;           // neighbour box id
  double dx_box;   // particle offset
  double dy_box;   // particle offset
  long k;          //particle index
  double Force_x,Force_y,force;    //Force between particles
  double dx,dy,dx2,dy2,dist,dist2,dist3;       // distance between particles
  
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

    //Loop through all particles k in box n
    k=Boxes[nbi][nbj]; //Start with the first
  
    //As long as there are particles in the box, iterate
    while(k!=-1){

      dx = Particles[j].x - Particles[k].x - dx_box; // positive if x_j>x_k
      dy = Particles[j].y - Particles[k].y - dy_box; // positive if y_j>y_k
      dx2=dx*dx;
      dy2=dy*dy;
      dist2=dx2+dy2;
    
      //If the particles are closer than Param.rmax, they interact
      if (dist2<Param.rmax2){
        dist=pow(dist2,0.5);
        force=Force_harmonic(dist, Param);
        Force_x = force*dx/dist;
        Force_y = force*dy/dist;
        forces[2*j] += Force_x;
        forces[2*j+1] += Force_y;
        forces[2*k] -= Force_x;
        forces[2*k+1] -= Force_y;
      }
      k=Neighbours[2*k+1]; // k is now the next particle in the list
    }
  }
}

void Loop_Force(particle* Particles, param Param, double* forces,long* Neighbours, 
                box*** NeighbouringBoxes, long** Boxes){
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
        Compute_force_same_box(j,Neighbours,forces,Param,Particles);

        // Loop through the neigboring box
        Compute_force_neighbours(j, bi, bj, NeighbouringBoxes, Particles, Param, forces, Neighbours, Boxes);
        
        //We have included all the interactions between j and other particles. Now iterate over j
        j=Neighbours[2*j+1];
      }
    // We are done with all particles in Box[i], iterate over the box
    }
  }
}


void Update_Particles(param Param,double* Displacements,particle* Particles, 
                      double* Integrated_Displacements, double* forces, long* Neighbours, 
                      box*** NeighbouringBoxes, long** Boxes,double _time){
  long n;
  int newbi,newbj;
  double v, cos_th, sin_th;
  
  // incorporate interactions with spatial hashing
  Loop_Force(Particles, Param, forces, Neighbours, NeighbouringBoxes, Boxes);
  
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

    // incorporate interactions
    Displacements[2*n]   += forces[2*n]*Param.dt;
    Displacements[2*n+1] += forces[2*n+1]*Param.dt;
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[2*n];
    Particles[n].y += Displacements[2*n+1];

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;

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
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[2*n];
    Particles[n].y += Displacements[2*n+1];

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;

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
  fprintf(output_param,"v0 is %lg\n", Param.v0);
  fprintf(output_param,"v1 is %lg\n", Param.v1);
  fprintf(output_param,"v_center1 is %lg\n", Param.v_center1);
  fprintf(output_param,"v_center2 is %lg\n", Param.v_center2);
  fprintf(output_param,"v_right is %lg\n", Param.v_right);
  fprintf(output_param,"P is %d\n", Param.P);
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
void Update_Profile(param Param, particle* Particles, long*** profile_rho){
  long i; // particle index
  int j,k;  // index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    profile_rho[0][j][k] += 1;
  }
}


// Store the profile of positions and polarizations
void Store_Profile(param Param, long*** profile_rho, double _time, FILE* output_profile_rho){
  int j,k; // profile index
  
  // columns: time, x, y, count
  for (j=0; j<Param.Nbinx; j++){
    for (k=0; k<Param.Nbiny; k++){
      fprintf(output_profile_rho, "%lg\t%lg\t%lg\t%ld\n", _time, ((double)j)*Param.binwidthx, ((double)k)*Param.binwidthy, profile_rho[0][j][k]);
      profile_rho[0][j][k]=0;
    }
  }
  fflush(output_profile_rho);
}
