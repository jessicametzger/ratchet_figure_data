


void Initialize_parameters(int argc, char* argv[],
  FILE** output_param, FILE** output_profile_rho, FILE** output_profile_m, FILE** output_pos, FILE** output_disp, FILE** output_traj,
  param* Param, double* _time, double* prev_percentage,
  particle** Particles,double** Displacements, double** Integrated_Displacements,
  long** profile_rho, long** profile_m,
  long long* seed,
  double* NextStoreProfile, double* NextUpdateProfile, double* NextStorePos, double* NextStoreDisp){
  
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
  strcat(command_base, "b1 "); argctarget ++;
  strcat(command_base, "b2 "); argctarget ++;
  strcat(command_base, "d "); argctarget ++;
  strcat(command_base, "P "); argctarget ++;
  strcat(command_base, "epsilon "); argctarget ++;
  strcat(command_base, "seed "); argctarget ++;
  strcat(command_base, "Nbin "); argctarget++;
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
  
  sprintf(name,"%s-profile_m",argv[i]);
  output_profile_m[0]=fopen(name,"w");

  sprintf(name,"%s-pos",argv[i]);
  output_pos[0]=fopen(name,"w");

  sprintf(name,"%s-disp",argv[i]);
  output_disp[0]=fopen(name,"w");
  
  i++;
  
  // Read in basic parameters
  Param[0].L                  = strtod(argv[i], NULL); i++;
  Param[0].N                  = (long) strtod(argv[i], NULL); i++;
  Param[0].dt                 = strtod(argv[i], NULL); i++;
  Param[0].final_time         = strtod(argv[i], NULL); i++;
  Param[0].alpha              = strtod(argv[i], NULL); i++;
  Param[0].v0                 = strtod(argv[i], NULL); i++;
  Param[0].a1                 = strtod(argv[i], NULL); i++;
  Param[0].a2                 = strtod(argv[i], NULL); i++;
  Param[0].b1                 = strtod(argv[i], NULL); i++;
  Param[0].b2                 = strtod(argv[i], NULL); i++;
  Param[0].d                  = strtod(argv[i], NULL); i++;
  Param[0].P                  = (int) strtod(argv[i], NULL); i++;
  Param[0].epsilon            = strtod(argv[i], NULL); i++;
  seed[0]                     = strtod(argv[i], NULL); i++;
  Param[0].Nbin               = strtod(argv[i], NULL); i++;
  Param[0].NstepProfile       = strtod(argv[i], NULL); i++;
  Param[0].StoreInterProfile  = strtod(argv[i], NULL); i++;
  Param[0].StoreInterPos      = strtod(argv[i], NULL); i++;
  Param[0].StoreInterDisp     = strtod(argv[i], NULL); i++;
  
  // define parameters that are functions of the inputs and/or known
  NextStoreProfile[0]         = Param[0].StoreInterProfile;
  NextUpdateProfile[0]        = NextStoreProfile[0] - (2/Param[0].alpha) * (Param[0].NstepProfile-1);
  NextStorePos[0]             = Param[0].StoreInterPos;
  NextStoreDisp[0]            = Param[0].StoreInterDisp;
  Param[0].coef1              = 2*M_PI*Param[0].P/Param[0].L;
  Param[0].coef2              = 4*M_PI*Param[0].P/Param[0].L;
  Param[0].shift1             = Param[0].d*2*M_PI;
  Param[0].shift1             = Param[0].d*4*M_PI;
  Param[0].v0_a1_half         = Param[0].v0*Param[0].a1/2.;
  Param[0].v0_a2_half         = Param[0].v0*Param[0].a2/2.;
  Param[0].potential_coef1    = Param[0].epsilon*Param[0].b1*(2*M_PI*Param[0].P/Param[0].L);
  Param[0].potential_coef2    = Param[0].epsilon*Param[0].b2*(4*M_PI*Param[0].P/Param[0].L);
  Param[0].alphadt            = Param[0].alpha*Param[0].dt;
  Param[0].binwidth           = (double) Param[0].L/((double) Param[0].Nbin);
  _time[0]                    = 0;
  prev_percentage[0]          = 0;
  
  // Must have enough time to make enough measurements between profile storages
  assert(Param[0].NstepProfile*(2/Param[0].alpha)<Param[0].StoreInterProfile);

  // Allocate space for simulation data
  Particles[0]                = (particle*) malloc(sizeof(particle)*Param[0].N);
  Displacements[0]            = (double*) calloc(Param[0].N,sizeof(double));
  Integrated_Displacements[0] = (double*) calloc(Param[0].N,sizeof(double));
  profile_rho[0]              = (long*) calloc(Param[0].Nbin,sizeof(long));
  profile_m[0]                = (long*) calloc(Param[0].Nbin,sizeof(long));
  
  // clear out profile
  memset(profile_rho[0],0,Param[0].Nbin*sizeof(long));
  memset(profile_m[0],0,Param[0].Nbin*sizeof(long));

// initialize particle locations and polarizations  
  for (i=0;i<Param[0].N;i++){
  
    // random initialize positions uniformly over interval
    Particles[0][i].x=Param[0].L*genrand64_real3();
    
    // randomly pick right or left motion
    Particles[0][i].theta = (genrand64_real2()<.5)?(-1):1;
    
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
  return Param.v0 + Param.v0_a1_half*sin(Param.coef1*x + Param.shift1) + Param.v0_a2_half*sin(Param.coef2*x + Param.shift2);
}

// Force_potential computes the force exerted on a particle by the repulsive potential epsilon/v.
// It returns epsilon*v'/v^2 = -U'(x).
double Force_potential(double x, param Param){
  return -Param.potential_coef1*cos(Param.coef1*x) - Param.potential_coef2*cos(Param.coef2*x);
}


void Update_Particles(param Param,double* Displacements,particle* Particles, double* Integrated_Displacements, double _time){
  long n;
  long k;
  double v;
  double minus_Uprime;
  
  for ( n=0 ; n<Param.N ; n++ ){
    // If the particle tumbles, it acquires a new direction
    if ( Particles[n].next_time < _time + Param.dt){
      Particles[n].theta = (Particles[n].theta>0)?(-1):1;
      Particles[n].next_time += -2*log(genrand64_real3())/Param.alpha;
    }
    
    // Compute self-propulsion velocity at this position
    v = vofx_function(Particles[n].x, Param);
    
    // Compute potential velocity at this postion
    minus_Uprime = Force_potential(Particles[n].x, Param);
    
    // Compute the displacement of each particle as x += v theta dt where theta= +/- 1
    Displacements[n] = Param.dt * Particles[n].theta * v + Param.dt * minus_Uprime;

  }
  
  for ( n=0 ; n<Param.N ; n++ ){
  
    // Once all displacement are computed, move the particles
    Particles[n].x += Displacements[n];

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.L)
      Particles[n].x-=Param.L;
    while(Particles[n].x<0)
      Particles[n].x+=Param.L;

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
  fprintf(output_param,"b1 is %lg\n", Param.b1);
  fprintf(output_param,"b2 is %lg\n", Param.b2);
  fprintf(output_param,"d is %lg\n", Param.d);
  fprintf(output_param,"P is %d\n", Param.P);
  fprintf(output_param,"epsilon is %lg\n", Param.epsilon);
  fprintf(output_param,"Nbin is %d\n", Param.Nbin);
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

// Update density and magnetization profiles (add one to bin for each particle currently in it)
void Update_Profile(param Param, particle* Particles, long* profile_rho, long* profile_m){
  long i; // particle index
  int j;  // index to add to. From 0 to Param.Nbin-1
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidth)+EPS);
    profile_rho[j] += 1;
    profile_m[j]   += Particles[i].theta;
  }
}


// Store the profile of positions and polarizations
void Store_Profile(param Param,particle* Particles, long** profile_rho, long** profile_m, double _time, FILE* output_profile_rho, FILE* output_profile_m){
  int j; // profile index
  
  for (j=0; j<Param.Nbin; j++){
    fprintf(output_profile_rho, "%lg\t%lg\t%ld\n", _time, ((double)j)*Param.binwidth, profile_rho[0][j]);
    fprintf(output_profile_m, "%lg\t%lg\t%ld\n", _time, ((double)j)*Param.binwidth, profile_m[0][j]);
  }
  fflush(output_profile_rho);
  fflush(output_profile_m);
      
  // clear out profile
  memset(profile_rho[0],0,Param.Nbin*sizeof(long));
  memset(profile_m[0],0,Param.Nbin*sizeof(long));
}