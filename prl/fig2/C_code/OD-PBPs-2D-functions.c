
void Initialize_parameters (int argc, char *argv[], param* Param, state* State, data* Data, particle** Particles, 
  long ***Boxes, long **Neighbours){

  long i, j, k, l;              // indices for allocating arrays
  double tmp;
  int argctarget = 0;           // Number of parameters that should be used
  char command_base[1000] = ""; // string that contains the desired format of the command line
  char name[200];               // string in which the file names are written

  /* Format the command line and count the arguments */
  argctarget = 0;
  strcat (command_base, "usage: ");
  strcat (command_base, argv[0]);                    argctarget ++;
  strcat (command_base, " ");
  strcat (command_base, "file ");                    argctarget ++;
  strcat (command_base, "Lx ");                      argctarget ++;
  strcat (command_base, "Ly ");                      argctarget ++;
  strcat (command_base, "N ");                       argctarget ++;
  strcat (command_base, "dt ");                      argctarget ++;
  strcat (command_base, "tf ");                      argctarget ++;
  strcat (command_base, "seed ");                    argctarget ++;
  strcat(command_base, "Nstep_T ");                  argctarget ++;
  strcat(command_base, "T_xs ");                     argctarget ++;
  strcat(command_base, "Ts ");                       argctarget ++;
  strcat (command_base, "a ");                       argctarget ++;
  strcat (command_base, "k ");                       argctarget ++;
  strcat (command_base, "StoreInterPos ");           argctarget ++;
  strcat(command_base, "Nbinx ");                    argctarget ++;
  strcat(command_base, "Nbiny ");                    argctarget ++;
  strcat(command_base, "NstepProf ");                argctarget ++;
  strcat(command_base, "StoreInterProf ");           argctarget ++;
  strcat(command_base, "UpdateInterProf ");          argctarget ++;
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
  Param[0].output_param = fopen (name, "w");

  sprintf (name, "%s-pos", argv[i]);
  Param[0].output_pos = fopen (name, "w");

  sprintf (name, "%s-Tprof", argv[i]);
  Param[0].output_Tprof = fopen (name, "w");
  
  sprintf(name,"%s-sigmaIKprof",argv[i]);
  Param[0].output_sigmaIKprof=fopen(name,"w");
  
  sprintf(name,"%s-EPR",argv[i]);
  Param[0].output_EPR=fopen(name,"w");

  printf ("created files\n");

  i++;

  // Read in basic parameters
  Param[0].Lx                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Ly                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].N                  = (long) strtod (argv[i], NULL)      ; i++ ;
  Param[0].dt                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].tf                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].seed               = (long long) strtod (argv[i], NULL) ; i++ ;
  Param[0].TParam.Nstep       = (int) strtod(argv[i], NULL)        ; i++ ;
  Param[0].Txs                = argv[i]                            ; i++ ;
  Param[0].Ts                 = argv[i]                            ; i++ ;
  Param[0].a                  = strtod (argv[i], NULL)             ; i++ ;
  Param[0].k                  = strtod (argv[i], NULL)             ; i++ ;
  Param[0].StoreInterPos      = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Nbinx              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Nbiny              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].NstepProf          = (int)strtod (argv[i], NULL)        ; i++ ;
  Param[0].StoreInterProf     = strtod (argv[i], NULL)             ; i++ ;
  Param[0].UpdateInterProf    = strtod (argv[i], NULL)             ; i++ ;
  

  // define parameters that are functions of the inputs and/or known

  Particles[0]             = (particle *)malloc (sizeof (particle) * Param[0].N);
  init_genrand64 (Param[0].seed);
  State[0]._time           = 0;
  State[0].prev_percentage = 0;
  State[0].NextStorePos    = Param[0].StoreInterPos;
  Param[0].UpdateInterProg = Param[0].tf/100.0; // update every 1%
  State[0].NextUpdateProg  = Param[0].UpdateInterProg;
  State[0].NextStoreProf   = Param[0].StoreInterProf;
  State[0].NextUpdateProf  = State[0].NextStoreProf - (Param[0].UpdateInterProf) * (Param[0].NstepProf - 1);
  assert((Param[0].NstepProf-1)*Param[0].UpdateInterProf <= Param[0].StoreInterProf);

  Param[0].sqrt2dt              = pow(2*Param[0].dt, 0.5);
  Param[0].TParam.L             = Param[0].Lx;
  Param[0].TParam.xs            = (double*) calloc(Param[0].TParam.Nstep, sizeof(double));
  Param[0].TParam.fs            = (double*) calloc(Param[0].TParam.Nstep, sizeof(double));
  split_vec(Param[0].TParam.Nstep, &Param[0].TParam.xs, Param[0].Txs, ',');
  split_vec(Param[0].TParam.Nstep, &Param[0].TParam.fs, Param[0].Ts, ',');
  compute_coefs_cub(&Param[0].TParam);

  
  Param[0].rmax            = Param[0].a;
  Param[0].rmax2           = Param[0].rmax*Param[0].rmax;
  Param[0].NxBox           = (int) (floor(Param[0].Lx/Param[0].rmax)+EPS);
  Param[0].NyBox           = (int) (floor(Param[0].Ly/Param[0].rmax)+EPS);

  // System width must be integer multiple of interaction length
  assert((double)Param[0].NxBox - Param[0].Lx/Param[0].rmax<EPS);
  assert((double)Param[0].NyBox - Param[0].Ly/Param[0].rmax<EPS);
  
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
  Param[0].NeighbouringBoxes = (box***) malloc(Param[0].NxBox*sizeof(box**));
  DefineNeighbouringBoxes(Param[0].NeighbouringBoxes,Param[0].NxBox, Param[0].NyBox, Param[0].Lx, Param[0].Ly);

  get_force_harmonic_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, 
                            &Param[0].c1, &Param[0].c2);



  Param[0].binwidthx = (double)Param[0].Lx / ((double)Param[0].Nbinx);
  Param[0].binwidthy = (double)Param[0].Ly / ((double)Param[0].Nbiny);

  Data[0].Tprof = calloc(Param[0].Nbinx, sizeof(double*));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].Tprof[i] = calloc(Param[0].Nbiny, sizeof(double));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].Tprof[i][j] = 0.0;
    }
  }

  Data[0].sigmaIKprof = calloc(Param[0].NxBox, sizeof(double**));
  for (i=0; i<Param[0].NxBox; i++){
    Data[0].sigmaIKprof[i] = calloc(Param[0].NyBox, sizeof(double*));
    for (j=0; j<Param[0].NyBox; j++){
      Data[0].sigmaIKprof[i][j] = calloc(4, sizeof(double));
      Data[0].sigmaIKprof[i][j][0] = 0.0; // sigmaxx
      Data[0].sigmaIKprof[i][j][1] = 0.0; // sigmaxy
      Data[0].sigmaIKprof[i][j][2] = 0.0; // sigmayy
    }
  }

  Data[0].EPR = calloc(Param[0].Nbinx, sizeof(double*));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].EPR[i] = calloc(Param[0].Nbiny, sizeof(double));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].EPR[i][j] = 0.0;
    }
  }


// initialize particle locations
  for (i = 0; i < Param[0].N; i++) {

      // random initialize positions uniformly over interval
      Particles[0][i].x = Param[0].Lx * genrand64_real2 ();
      Particles[0][i].y = Param[0].Ly * genrand64_real2 ();

      Particles[0][i].dx = 0.0;
      Particles[0][i].dy = 0.0;

      // Add particle to its box and link to its neighbors
      Particles[0][i].bi = (int)(floor (Particles[0][i].x / Param[0].rmax) + EPS);
      Particles[0][i].bj = (int)(floor (Particles[0][i].y / Param[0].rmax) + EPS);
      AddinBox (i, Particles[0][i].bi, Particles[0][i].bj, Boxes[0], Neighbours[0]);

      Particles[0][i].ifx = 0.0;
      Particles[0][i].ify = 0.0;
      Particles[0][i].ifxs = (double*) calloc(Param[0].N, sizeof(double));
      Particles[0][i].ifys = (double*) calloc(Param[0].N, sizeof(double));
      for (j=0; j<Param[0].N; j++) {
        Particles[0][i].ifxs[j] = 0.0;
        Particles[0][i].ifys[j] = 0.0;
      }
    }
}

double Tofx(fparam_cub TParam, double x, double y){
  return fofx_cub(TParam, x);
}

void Update_Particles(param Param,particle* Particles, state State, long* Neighbours, long** Boxes){
  long n;
  int newbi,newbj;
  double T,sqrt2Tdt;

  // reset displacements and forces to 0
  for (n=0; n<Param.N; n++){
    Particles[n].dx = 0.0;
    Particles[n].dy = 0.0;
    Particles[n].ifx = 0.0;
    Particles[n].ify = 0.0;
    memset(Particles[n].ifxs, 0.0, Param.N);
    memset(Particles[n].ifys, 0.0, Param.N);
  }
  
  // incorporate interactions
  Loop_Force(Param, Particles, Neighbours, Boxes);

  for (n=0 ; n<Param.N ; n++ ){

    // Compute the displacement of each particle

    Particles[n].T = Tofx(Param.TParam, Particles[n].x, Particles[n].y);
    sqrt2Tdt = pow(Particles[n].T, 0.5) * Param.sqrt2dt;

    Particles[n].dx += sqrt2Tdt * gasdev();
    Particles[n].dy += sqrt2Tdt * gasdev();

    // incorporate interactions
    Particles[n].dx += Particles[n].ifx;
    Particles[n].dy += Particles[n].ify;
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
    // Once all displacement are computed, move the particles
    Particles[n].x += Particles[n].dx;
    Particles[n].y += Particles[n].dy;

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;

    // Update box membership
    newbi = (int) (floor(Particles[n].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[n].y/Param.rmax) + EPS);

    if(Particles[n].bi!=newbi || Particles[n].bj!=newbj){
      RemovefromBox(n,Particles[n].bi,Particles[n].bj,Boxes,Neighbours);
      AddinBox(n,newbi,newbj,Boxes,Neighbours);
      Particles[n].bi=newbi;
      Particles[n].bj=newbj;
    }
  }
}


void Store_Parameters (int argc, char *argv[], param Param) {
  long i;

  /* Store parameters */
  for (i = 0; i < argc; i++) fprintf (Param.output_param, "%s ", argv[i]);

  fprintf (Param.output_param, "\n");
  fprintf (Param.output_param, "Lx is %lg\n",              Param .Lx);
  fprintf (Param.output_param, "Ly is %lg\n",              Param .Ly);
  fprintf (Param.output_param, "N is %ld\n",               Param .N);
  fprintf (Param.output_param, "dt is %lg\n",              Param .dt);
  fprintf (Param.output_param, "tf is %lg\n",              Param .tf);
  fprintf (Param.output_param, "seed is %lld\n",           Param .seed);
  fprintf( Param .output_param, "Nstep_T is %d\n",         Param .TParam .Nstep);
  fprintf( Param .output_param, "Txs is %s\n",             Param .Txs);
  fprintf( Param .output_param, "Ts is %s\n",              Param .Ts);
  fprintf (Param.output_param, "a is %lg\n",               Param .a);
  fprintf (Param.output_param, "k is %lg\n",               Param .k);
  fprintf (Param.output_param, "StoreInterPos is %lg\n",   Param .StoreInterPos);
  fprintf (Param.output_param, "Nbinx is %d\n",            Param .Nbinx);
  fprintf (Param.output_param, "Nbiny is %d\n",            Param .Nbiny);
  fprintf (Param.output_param, "NstepProf is %d\n",        Param .NstepProf);
  fprintf (Param.output_param, "StoreInterProf is %lg\n",  Param .StoreInterProf);
  fprintf (Param.output_param, "UpdateInterProf is %lg\n", Param .UpdateInterProf);
  fflush (Param.output_param);
}


// Store the identities, positions of all particles, with the time
void Store_Pos(param Param,particle* Particles,state State){

  long n;
  fprintf(Param.output_pos, "%.4f",State._time);
  for (n=0;n<Param.N;n++){
    fprintf(Param.output_pos,"\t%lg",Particles[n].x);
  }
  fprintf(Param.output_pos,"\n%.4f",State._time);
  for (n=0;n<Param.N;n++){
    fprintf(Param.output_pos,"\t%lg",Particles[n].y);
  }
  fprintf(Param.output_pos,"\n");
  fflush(Param.output_pos);
}


void Update_T(param Param, particle* Particles, data* Data){
  long i;
  int j,k;
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    Data[0].Tprof[j][k] += Particles[i].T;
  }
}
void Store_T(param Param,particle* Particles, state State, data* Data){
  int n,j,k; // profile index
  
  for (j=0; j<Param.Nbinx; j++){
    fprintf(Param.output_Tprof, "%.4f",State._time);
    for (k=0; k<Param.Nbiny; k++){
      fprintf(Param.output_Tprof, "\t%lg", Data[0].Tprof[j][k]);
      Data[0].Tprof[j][k]=0.0;
    }
    fprintf(Param.output_Tprof, "\n");
  }
  fflush(Param.output_Tprof);
}


void Store_SigmaIK(param Param,particle* Particles, state State, data* Data) {
  int n,j,k; // profile index
  
  for (n=0; n<3; n++){
    for (j=0; j<Param.NxBox; j++){
      fprintf(Param.output_sigmaIKprof, "%.4f",State._time);
      for (k=0; k<Param.NyBox; k++){
        fprintf(Param.output_sigmaIKprof, "\t%lg", Data[0].sigmaIKprof[j][k][n]);
        Data[0].sigmaIKprof[j][k][n]=0.0;
      }
      fprintf(Param.output_sigmaIKprof, "\n");
    }
  }
  fflush(Param.output_sigmaIKprof);
}




// calculate EPR -- requires Stratonovich average, so moving particles back by dx/2
void Update_EPR(param Param, particle* Particles, data* Data, long* Neighbours, long** Boxes){
  long i;  // particle index
  int j,k; // profile index
  int newbi,newbj; // box membership

  // move to position at middle of timestep
  for (i=0; i<Param.N; i++){
    Particles[i].x -= Particles[i].dx * 0.5;
    Particles[i].y -= Particles[i].dy * 0.5;

    Particles[i].ifx = 0.0;
    Particles[i].ify = 0.0;

    // Take care of periodic boundary conditions
    while(Particles[i].x>Param.Lx) Particles[i].x-=Param.Lx;
    while(Particles[i].x<0) Particles[i].x+=Param.Lx;
    while(Particles[i].y>Param.Ly) Particles[i].y-=Param.Ly;
    while(Particles[i].y<0) Particles[i].y+=Param.Ly;


    // Update box membership
    newbi = (int) (floor(Particles[i].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[i].y/Param.rmax) + EPS);

    if(Particles[i].bi!=newbi || Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Boxes,Neighbours);
      AddinBox(i,newbi,newbj,Boxes,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
  }

  // calculate interaction force at middle of timestep
  Loop_Force(Param, Particles, Neighbours, Boxes);

  // calculate EPR
  for (i=0; i<Param.N; i++){

    // calculate temperature at intermediate location
    Particles[i].T = Tofx(Param.TParam, Particles[i].x, Particles[i].y);

    j = (int) (floor(Particles[i].x/Param.binwidthx) + EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy) + EPS);
    Data[0].EPR[j][k] += (Particles[i].dx * Particles[i].ifx + Particles[i].dy * Particles[i].ify) / (Particles[i].T * Param.dt);
  }

  // move back to original position 
  for (i=0; i<Param.N; i++){
    Particles[i].x += Particles[i].dx * 0.5;
    Particles[i].y += Particles[i].dy * 0.5;

    // Take care of periodic boundary conditions
    while(Particles[i].x>Param.Lx) Particles[i].x-=Param.Lx;
    while(Particles[i].x<0) Particles[i].x+=Param.Lx;
    while(Particles[i].y>Param.Ly) Particles[i].y-=Param.Ly;
    while(Particles[i].y<0) Particles[i].y+=Param.Ly;

    // Update box membership
    newbi = (int) (floor(Particles[i].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[i].y/Param.rmax) + EPS);

    if(Particles[i].bi!=newbi || Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Boxes,Neighbours);
      AddinBox(i,newbi,newbj,Boxes,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
  }
}

void Store_EPR(param Param, particle* Particles, state State, data* Data){
  int i,j; // profile index

  for (i=0; i<Param.Nbinx; i++){
    fprintf(Param.output_EPR, "%.4f", State._time);
    for (j=0; j<Param.Nbiny; j++){
      fprintf(Param.output_EPR, "\t%lg", Data[0].EPR[i][j]);
      Data[0].EPR[i][j] = 0.0;
    }
    fprintf(Param.output_EPR,"\n");
  }
  fflush(Param.output_EPR);
}



void Update_Data(param Param, particle* Particles, state* State, data* Data, long* Neighbours, long** Boxes){
  if (State[0]._time > State[0].NextUpdateProf-EPS){
    Update_T(Param, Particles, Data);
    Update_SigmaIK(Param, Particles, Data);
    Update_EPR(Param, Particles, Data, Neighbours, Boxes);
    State[0].NextUpdateProf += Param.UpdateInterProf;
  }
}

void Store_Data(param Param, particle* Particles, state* State, data* Data){
  if (State[0]._time > State[0].NextStoreProf-EPS){
    Store_T(Param, Particles, State[0], Data);
    Store_SigmaIK(Param, Particles, State[0], Data);
    Store_EPR(Param, Particles, State[0], Data);
    State[0].NextStoreProf += Param.StoreInterProf;
  }
}


/* Print percentage of completed simulation */
void Update_Progress(param Param, state* State){
  printf("\rIn progress [%d %%]", (int)(State[0].prev_percentage + 1.01));
  fflush(stdout);
  State[0].prev_percentage += 1.0;
}