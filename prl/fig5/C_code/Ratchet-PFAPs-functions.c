
// Add particle i into box bi and set its new neighbors. Make new particle the 1st in box
void AddinBox(long i, int bi, long* Boxes, long* Neighbors){
  long k;
  
  // Save old 1st particle in box
  // Or, if there are no particles in bi, k=-1
  k = Boxes[bi];
  
  // Make new particle new 1st particle in box
  Boxes[bi] = i;
  Neighbors[2*i]=-1;
  
  // Update linked list so k (original 1st) follows i.
  Neighbors[2*i+1]=k;
  if (k!=-1) {
    Neighbors[2*k] = i;
  }
}

// Remove particle i from box bi and link original neighbors.
// bi is index of original box (box it's being removed from)
void RemovefromBox(long i, int bi, long* Boxes, long* Neighbors){
  long next;
  long prev;
  
  // Store next particle
  next = Neighbors[2*i+1];
  
  // If particle is 1st in box
  if (Boxes[bi]==i){
    // Make the 2nd particle in box the new 1st
    Boxes[bi]=next;
    
    // Make the 2nd particle (now, first) have no one before it
    if (next!=-1){
      Neighbors[2*next]=-1;
    }
  }
  // If the particle isn't 1st in box
  else {
    // Store previous particle
    prev=Neighbors[2*i];
    
    // Link the previous to the enext
    Neighbors[2*prev+1]=next;
    
    // Make the next particle's previous our original previous
    if (next!=-1){
      Neighbors[2*next]=prev;
    }
  }
}

// initialize empty boxes
void ConstructBoxes(param Param, long* Boxes){
  int i; // box index
  for (i=0; i<Param.Nbox; i++){
    Boxes[i] = -1;
  }
}

// construct the neighbor array. start it as no neighbors.
void ConstructNeighbors(param Param, long* Neighbors){
  long i; // neighbor index
  for (i=0; i<Param.N; i++){
    Neighbors[i] = -1;
  }
}

// construct list of next box
void ConstructNextBoxes(param Param, box* NextBoxes){
  int i;
  for (i=0; i<Param.Nbox; i++){
    NextBoxes[i].index = (i+1);
    NextBoxes[i].epsilonx = 0;
  }
#ifdef PBC
  NextBoxes[Param.Nbox-1].index = 0;
  NextBoxes[Param.Nbox-1].epsilonx = Param.L;
#endif 
}


void Initialize_parameters(int argc, char* argv[],
  FILE** output_param, FILE** output_pos, FILE** output_disp,
  param* Param, double* _time, double* prev_percentage,
  particle** Particles,double** Displacements, double** Integrated_Displacements, double** forces, long** Boxes, long** Neighbors, box** NextBoxes,
  long long* seed,double* NextStorePos, double* NextStoreDisp){
  
  long i;
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
  strcat(command_base, "L "); argctarget ++;
  strcat(command_base, "N "); argctarget ++;
  strcat(command_base, "dt "); argctarget ++;
  strcat(command_base, "final_time "); argctarget ++;
  strcat(command_base, "alpha "); argctarget ++;
  strcat(command_base, "v0 "); argctarget ++;
  strcat(command_base, "a1 "); argctarget ++;
  strcat(command_base, "a2 "); argctarget ++;
  strcat(command_base, "P "); argctarget ++;
  strcat(command_base, "sigma "); argctarget ++;
  strcat(command_base, "amp "); argctarget ++;
  strcat(command_base, "seed "); argctarget ++;
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

  sprintf(name,"%s-pos",argv[i]);
  output_pos[0]=fopen(name,"w");

  sprintf(name,"%s-disp",argv[i]);
  output_disp[0]=fopen(name,"w");
  
  i++;
  
  // Read in basic parameters
  Param[0].L                 = strtod(argv[i], NULL); i++;
  Param[0].N                 = (long) strtod(argv[i], NULL); i++;
  Param[0].dt                = strtod(argv[i], NULL); i++;
  Param[0].final_time        = strtod(argv[i], NULL); i++;
  Param[0].alpha             = strtod(argv[i], NULL); i++;
  Param[0].v0                = strtod(argv[i], NULL); i++;
  Param[0].a1                = strtod(argv[i], NULL); i++;
  Param[0].a2                = strtod(argv[i], NULL); i++;
  Param[0].P                 = (int) strtod(argv[i], NULL); i++;
  Param[0].sigma             = strtod(argv[i], NULL); i++;
  Param[0].amp               = strtod(argv[i], NULL); i++;
  seed[0]                    = strtod(argv[i], NULL); i++;
  Param[0].StoreInterPos     = strtod(argv[i], NULL); i++;
  Param[0].StoreInterDisp    = strtod(argv[i], NULL); i++;
  
  // define parameters that are functions of the inputs and/or known
  NextStorePos[0]      = Param[0].StoreInterPos;
  NextStoreDisp[0]     = Param[0].StoreInterDisp;
  Param[0].coef1       = 2*M_PI*Param[0].P/Param[0].L;
  Param[0].coef2       = 4*M_PI*Param[0].P/Param[0].L;
  Param[0].v0_a1_half  = Param[0].v0*Param[0].a1/2.;
  Param[0].v0_a2_half  = Param[0].v0*Param[0].a2/2.;
  Param[0].rmax        = Param[0].sigma;
  Param[0].rmax2       = Param[0].rmax*Param[0].rmax;
  Param[0].max_amp     = Param[0].amp/Param[0].sigma;
  Param[0].alphadt     = Param[0].alpha*Param[0].dt;
  Param[0].rbox        = Param[0].rmax;
  Param[0].Nbox        = (int) (Param[0].L / Param[0].rbox + EPS);
  _time[0]             = 0;
  prev_percentage[0]   = 0;
  
  init_genrand64(seed[0]);
  
  // System width must be integer multiple of interaction length
  assert((double)Param[0].Nbox - Param[0].L/Param[0].rbox<EPS);

  // Allocate space for simulation data
  Particles[0]                = (particle*) malloc(sizeof(particle)*Param[0].N);
  Displacements[0]            = (double*) calloc(Param[0].N,sizeof(double));
  Integrated_Displacements[0] = (double*) calloc(Param[0].N,sizeof(double));
  forces[0]                   = (double*) calloc(Param[0].N,sizeof(double));
  
  // Construct list of 1st particle in each box (initialized as empty)
  Boxes[0] = (long*) calloc(Param[0].Nbox, sizeof(long));
  ConstructBoxes(Param[0], Boxes[0]);
  
  // Construct empty list of neighbors
  Neighbors[0] = (long*) calloc((Param[0].N+1)*2,sizeof(long));
  ConstructNeighbors(Param[0], Neighbors[0]);
  
  // Construct list of next boxes
  NextBoxes[0] = (box*) malloc(Param[0].Nbox*sizeof(box));
  ConstructNextBoxes(Param[0], NextBoxes[0]);

// initialize particle locations and polarizations  
  for (i=0;i<Param[0].N;i++){
  
    // random initialize positions uniformly over interval
    Particles[0][i].x=Param[0].L*genrand64_real3();
    
    // randomly pick right or left motion
    Particles[0][i].theta = (genrand64_real2()<.5)?(-1):1;
    
    // Add particle to its box and link to its neighbors
    Particles[0][i].bi= (int) (floor(Particles[0][i].x/Param[0].rbox) + EPS);
    AddinBox(i,Particles[0][i].bi,Boxes[0],Neighbors[0]);
    
    // initialize next flip time from exponential distribution
    Particles[0][i].next_time = -2*log(genrand64_real3())/Param[0].alpha;
  }
}

/*
  There are two choices for space-dependent speed: a function that one
  calls or an array that one reads. We should try both.
*/

double vofx_function(double x, param Param){
  /*
    Make a periodic velocity field made of two Fourier modes.
  */

  return Param.v0 + Param.v0_a1_half*sin(Param.coef1*x) + Param.v0_a2_half*sin(Param.coef2*x);
}

// Force_harmonic computes the repulsive harmonic force exerted by the particle k onto particle j
// if dx>0, the force has to be positive
double Force_harmonic(double dx, param Param){
  return Param.max_amp*dx;
}


//This function computes the interactions between particle J and all
//the subsequent particles in the same box
void Compute_force_same_box(long j,long* Neighbors,double* forces, param Param, particle* Particles){
  long k; // index of neighbors
  double Force; // Force exerted between two particles
  double dx; // displacement between particles
  double dist2; // distance between particles squared
  
  // Loop through all the following neighbors of j
  k=Neighbors[2*j+1]; // start with the 1st neighbor
  
  // As long as there are neighbors, iterate
  while(k!=-1){
  
    dx=Particles[j].x-Particles[k].x;// positive if x_j>x_k
    Force=Force_harmonic(dx,Param);
    forces[j] += Force;
    forces[k] -= Force;
    k=Neighbors[2*k+1]; // k is now the next particle in the list
  }
}

// This functions compute the interactions between particle j and all
// the particles in the box situated on its right, using periodic
// boundary conditions
void Compute_force_neighbors(long j,int bi,box* NextBoxes,particle* Particles,param Param,double* forces, long* Neighbors, long* Boxes){
  int n;           //box name
  double dx_box;   // particle offset
  long k;          //particle index
  double dist2;    // Squared distance between particles
  double Force;    //Force between particles
  double dx;       // distance between particles
  
  // In 1d, there is only one neighboring box to the box i in which j is
  n=NextBoxes[bi].index;        // name of the box
  dx_box=NextBoxes[bi].epsilonx;// offset to be added to the particles in
			         // box bi to take into account the periodic
			         // boundary conditions

  //Loop through all particles k in box n
  k=Boxes[n]; //Start with the first
  
  //As long as there are particles in the box, iterate
  while(k!=-1){
    dx=Particles[j].x-Particles[k].x-dx_box;// negative always
    dist2=dx*dx;
    //If the particles are closer than Param.rmax, they interact
    if (dist2<Param.rmax2){
      Force=Force_harmonic(dx,Param);
      forces[j] += Force;
      forces[k] -= Force;
    }
    k=Neighbors[2*k+1]; // k is now the next particle in the list
  }
}

void Loop_Force(particle* Particles, param Param, double* forces,long* Neighbors, box* NextBoxes, long* Boxes){
  int i; // Box index
  long j; // Particle index

  //Initialize the array of forces to zero
  memset(forces,0,Param.N*sizeof(double));
  
  /* 
     To compute the force on each particle, we loop through all the boxes and compute:
     - the interactions between particles inside the box
     - the interactions between particles inside the box and inside a neighboring box
  */

  // Loop through all boxes
  for(i=0; i<Param.Nbox; i++){
  
    //Compute the force inside the box

    //Loop through all the particles 'j' in the box.
    j=Boxes[i]; //Start with j being the 1st particle

    //As long as j is not -1, compute its interactions with all the other particles
    while(j!=-1){
    
      // Loop first with the particles in the same box
      Compute_force_same_box(j,Neighbors,forces,Param,Particles);

      // Loop through the neigboring box
      Compute_force_neighbors(j, i, NextBoxes, Particles, Param, forces, Neighbors, Boxes);
      
      //We have included all the interactions between j and other particles. Now iterate over j
      j=Neighbors[2*j+1];
    }
    // We are done with all particles in Box[i], iterate over the box
  }
}


void Update_Particles(param Param,double* Displacements,particle* Particles, 
    double* Integrated_Displacements, double* forces, long* Neighbors, box* NextBoxes, 
    long* Boxes, double _time){
  long n;
  int newbi;
  
  
  // incorporate interactions with spatial hashing
  Loop_Force(Particles, Param, forces, Neighbors, NextBoxes, Boxes);
  
  for ( n=0 ; n<Param.N ; n++ ){
    // If the particle tumbles, it acquires a new direction
    if ( Particles[n].next_time < _time + Param.dt){
      Particles[n].theta = (Particles[n].theta>0)?(-1):1;
      Particles[n].next_time += -2*log(genrand64_real3())/Param.alpha;
    }
    // Compute the displacement of each particle as x += v theta dt where theta= +/- 1
    Displacements[n] = Param.dt * Particles[n].theta * vofx_function(Particles[n].x, Param);

    // incorporate interactions
    Displacements[n] += forces[n]*Param.dt;
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[n];

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.L)
      Particles[n].x-=Param.L;
    while(Particles[n].x<0)
      Particles[n].x+=Param.L;

    // Update box membership
    newbi = (int) (floor(Particles[n].x/Param.rbox) + EPS);

    if(Particles[n].bi!=newbi){
      RemovefromBox(n,Particles[n].bi,Boxes,Neighbors);
      AddinBox(n,newbi,Boxes,Neighbors);
      Particles[n].bi=newbi;
    }

    // Update the cumulated displacements
    Integrated_Displacements[n] += Displacements[n];
  }
}


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
  
  fprintf(output_param,"\n");
  
  fprintf(output_param,"L is %lg\n", Param.L);
  fprintf(output_param,"N is %ld\n", Param.N);
  fprintf(output_param,"dt is %lg\n", Param.dt);
  fprintf(output_param,"final_time is %lg\n", Param.final_time);
  fprintf(output_param,"alpha is %lg\n", Param.alpha);
  fprintf(output_param,"v0 is %lg\n", Param.v0);
  fprintf(output_param,"a1 is %lg\n", Param.a1);
  fprintf(output_param,"a2 is %lg\n", Param.a2);
  fprintf(output_param,"P is %d\n", Param.P);
  fprintf(output_param,"sigma is %lg\n", Param.sigma);
  fprintf(output_param,"amp is %lg\n", Param.amp);
  fprintf(output_param,"StoreInterPos is %lg\n", Param.StoreInterPos);
  fprintf(output_param,"StoreInterDisp is %lg\n", Param.StoreInterDisp);
  fflush(output_param);
}

// Store the identities, positions, and orientations of all particles, with the time
void Store_Pos(param Param,particle* Particles,double _time,FILE* output_pos){

  long n;
  for (n=0;n<Param.N;n++){
    fprintf(output_pos,"%lg\t%ld\t%lg\t%lg\n",_time,n,Particles[n].x,Particles[n].theta);
  }
  fflush(output_pos);
}

// Store the identities and integrated displacements of all particles, with the time
void Store_Disp(param Param, double* Integrated_Displacements, double _time,FILE* output_disp){

  long n;
  for (n=0;n<Param.N;n++){
    fprintf(output_disp,"%lg\t%ld\t%lg\n",_time,n,Integrated_Displacements[n]);
  }
  fflush(output_disp);
}