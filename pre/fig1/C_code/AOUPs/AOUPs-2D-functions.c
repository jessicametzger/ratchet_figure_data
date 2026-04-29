
void Initialize_parameters (int argc, char *argv[], param* Param, state* State, data* Data, particle** Particles
#ifndef NI
                       , long ***Boxes, long **Neighbours
#endif
                       ){

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
  strcat (command_base, "tau ");                     argctarget ++;
  strcat(command_base, "Nstep_D ");                  argctarget ++;
  strcat(command_base, "Dxs ");                      argctarget ++;
  strcat(command_base, "Ds ");                       argctarget ++;
#ifndef NI
  strcat (command_base, "a ");                       argctarget ++;
  strcat (command_base, "k ");                       argctarget ++;
#endif
  strcat (command_base, "StoreInterPos ");           argctarget ++;
  strcat(command_base, "StoreInterDisp ");           argctarget ++;
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

  sprintf (name, "%s-disp", argv[i]);
  Param[0].output_disp = fopen (name, "w");

  sprintf (name, "%s-prof", argv[i]);
  Param[0].output_prof = fopen (name, "w");

  sprintf (name, "%s-Jprof", argv[i]);
  Param[0].output_Jprof = fopen (name, "w");

  printf ("created files\n");
  fflush(stdout);

  i++;

  // Read in basic parameters
  Param[0].Lx                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Ly                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].N                  = (long) strtod (argv[i], NULL)      ; i++ ;
  Param[0].dt                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].tf                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].seed               = (long long) strtod (argv[i], NULL) ; i++ ;
  Param[0].tau                = strtod (argv[i], NULL)             ; i++ ;
  Param[0].DParam.Nstep       = (int) strtod(argv[i], NULL)        ; i++ ;
  Param[0].Dxs                = argv[i]                            ; i++ ;
  Param[0].Ds                 = argv[i]                            ; i++ ;
#ifndef NI
  Param[0].a                  = strtod (argv[i], NULL)             ; i++ ;
  Param[0].k                  = strtod (argv[i], NULL)             ; i++ ;
#endif
  Param[0].StoreInterPos      = strtod (argv[i], NULL)             ; i++ ;
  Param[0].StoreInterDisp     = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Nbinx              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Nbiny              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].NstepProf          = (int)strtod (argv[i], NULL)        ; i++ ;
  Param[0].StoreInterProf     = strtod (argv[i], NULL)             ; i++ ;
  Param[0].UpdateInterProf    = strtod (argv[i], NULL)             ; i++ ;
  

  // define parameters that are functions of the inputs and/or known

  Particles[0]             = (particle *)malloc (sizeof (particle) * Param[0].N);
  Param[0].dt_tau          = Param[0].dt/Param[0].tau;
  init_genrand64 (Param[0].seed);
  State[0]._time           = 0;
  State[0].prev_percentage = 0;
  State[0].NextStorePos    = Param[0].StoreInterPos;
  Param[0].UpdateInterProg = Param[0].tf/100.0; // update every 1%
  State[0].NextUpdateProg  = Param[0].UpdateInterProg;
  State[0].NextStoreDisp   = Param[0].StoreInterDisp;

  State[0].NextStoreProf   = Param[0].StoreInterProf;
  State[0].NextUpdateProf  = State[0].NextStoreProf - (Param[0].UpdateInterProf) * (Param[0].NstepProf - 1);
  assert((Param[0].NstepProf-1)*Param[0].UpdateInterProf <= Param[0].StoreInterProf);

  Param[0].binwidthx = (double)Param[0].Lx / ((double)Param[0].Nbinx);
  Param[0].binwidthy = (double)Param[0].Ly / ((double)Param[0].Nbiny);


  Param[0].sqrt2dt_tau = pow(2 * Param[0].dt, 0.5) / Param[0].tau;

  Param[0].DParam.L             = Param[0].Lx;
  Param[0].DParam.xs            = (double*) calloc(Param[0].DParam.Nstep, sizeof(double));
  Param[0].DParam.fs            = (double*) calloc(Param[0].DParam.Nstep, sizeof(double));
  split_vec(Param[0].DParam.Nstep, &Param[0].DParam.xs, Param[0].Dxs, ',');
  split_vec(Param[0].DParam.Nstep, &Param[0].DParam.fs, Param[0].Ds, ',');

  Param[0].DParam.cs            = (double*) calloc(Param[0].DParam.Nstep, sizeof(double));
  Param[0].DParam.ds            = (double*) calloc(Param[0].DParam.Nstep, sizeof(double));
  for (i=0; i<Param[0].DParam.Nstep-1; i++){
    tmp                         = Param[0].DParam.xs[i+1] - Param[0].DParam.xs[i];
    Param[0].DParam.cs[i]       = -3 * (Param[0].DParam.fs[i] - Param[0].DParam.fs[i+1]) / pow(tmp,2);
    Param[0].DParam.ds[i]       = 2 * (Param[0].DParam.fs[i] - Param[0].DParam.fs[i+1]) / pow(tmp,3);
  }
  tmp                           = Param[0].DParam.xs[0] + Param[0].Lx - Param[0].DParam.xs[i];
  Param[0].DParam.cs[i]         = -3 * (Param[0].DParam.fs[i] - Param[0].DParam.fs[0]) / pow(tmp,2);
  Param[0].DParam.ds[i]         = 2 * (Param[0].DParam.fs[i] - Param[0].DParam.fs[0]) / pow(tmp,3);


#ifndef NI
  
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
#endif


  Data[0].prof = calloc (Param[0].Nbinx, sizeof (long **));
  for (i = 0; i < Param[0].Nbinx; i++) {
    Data[0].prof[i] = calloc (Param[0].Nbiny, sizeof (long *));
    for (j = 0; j < Param[0].Nbiny; j++) {
      Data[0].prof[i][j] = 0.0;
    }
  }

  Data[0].Jprof = calloc(Param[0].Nbinx, sizeof(double**));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].Jprof[i] = calloc(Param[0].Nbiny, sizeof(double*));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].Jprof[i][j] = calloc(4, sizeof(double));
      Data[0].Jprof[i][j][0] = 0.0; // fx
      Data[0].Jprof[i][j][1] = 0.0; // fy
      Data[0].Jprof[i][j][2] = 0.0; // sqrt(2*D) eta_x
      Data[0].Jprof[i][j][3] = 0.0; // sqrt(2*D) eta_y
    }
  }

// initialize particle locations
  for (i = 0; i < Param[0].N; i++) {

      // random initialize positions uniformly over interval
      Particles[0][i].x = Param[0].Lx * genrand64_real2 ();
      Particles[0][i].y = Param[0].Ly * genrand64_real2 ();
      Particles[0][i].vx = 0.0;
      Particles[0][i].vy = 0.0;
      Particles[0][i].dx = 0.0;
      Particles[0][i].dy = 0.0;
      Particles[0][i].dvx = 0.0;
      Particles[0][i].dvy = 0.0;

      Particles[0][i].idx = 0.0;
      Particles[0][i].idy = 0.0;

#ifndef NI
      // Add particle to its box and link to its neighbors
      Particles[0][i].bi = (int)(floor (Particles[0][i].x / Param[0].rmax) + EPS);
      Particles[0][i].bj = (int)(floor (Particles[0][i].y / Param[0].rmax) + EPS);
      AddinBox (i, Particles[0][i].bi, Particles[0][i].bj, Boxes[0], Neighbours[0]);

      Particles[0][i].ifx = 0.0;
      Particles[0][i].ify = 0.0;
#endif
    }
}

double Dofx(fparam_cub DParam, double x, double y){
  return fofx_cub(DParam, x);
}


void Update_Particles(param Param,particle* Particles, state State
#ifndef NI
                      , long* Neighbours, long** Boxes
#endif
                      ){
  long n;
  int newbi,newbj;
  double D,sqrt2Ddt_tau;

  // reset displacements and forces to 0
  for (n=0; n<Param.N; n++){
    Particles[n].dx = 0.0;
    Particles[n].dy = 0.0;
    Particles[n].dvx = 0.0;
    Particles[n].dvy = 0.0;
#ifndef NI 
    Particles[n].ifx = 0.0;
    Particles[n].ify = 0.0;
#endif
  }
  
  // incorporate interactions
#ifndef NI
  Loop_Force(Param, Particles, Neighbours, Boxes);
#endif

  for (n=0 ; n<Param.N ; n++ ){

    // Compute the displacement of each particle
    Particles[n].D = Dofx(Param.DParam, Particles[n].x, Particles[n].y);
    sqrt2Ddt_tau = pow(Particles[n].D, 0.5) * Param.sqrt2dt_tau;

    Particles[n].dvx += sqrt2Ddt_tau * gasdev();
    Particles[n].dvy += sqrt2Ddt_tau * gasdev();
    Particles[n].dvx -= Param.dt_tau * Particles[n].vx;
    Particles[n].dvy -= Param.dt_tau * Particles[n].vy;

    Particles[n].dx = Particles[n].vx * Param.dt;
    Particles[n].dy = Particles[n].vy * Param.dt;

#ifndef NI
    // incorporate interactions
    Particles[n].dx += Particles[n].ifx;
    Particles[n].dy += Particles[n].ify;
#endif
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
    // Once all displacement are computed, move the particles
    Particles[n].x  += Particles[n].dx;
    Particles[n].y  += Particles[n].dy;
    Particles[n].vx += Particles[n].dvx;
    Particles[n].vy += Particles[n].dvy;

    // Take care of periodic boundary conditions
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;

#ifndef NI
    // Update box membership
    newbi = (int) (floor(Particles[n].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[n].y/Param.rmax) + EPS);

    if(Particles[n].bi!=newbi || Particles[n].bj!=newbj){
      RemovefromBox(n,Particles[n].bi,Particles[n].bj,Boxes,Neighbours);
      AddinBox(n,newbi,newbj,Boxes,Neighbours);
      Particles[n].bi=newbi;
      Particles[n].bj=newbj;
    }
#endif

    // Update the cumulated displacements
    Particles[n].idx += Particles[n].dx;
    Particles[n].idy += Particles[n].dy;
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
  fprintf (Param.output_param, "tau is %lg\n",             Param .tau);
  fprintf( Param .output_param, "Nstep_D is %d\n",         Param .DParam .Nstep);
  fprintf( Param .output_param, "Dxs is %s\n",             Param .Dxs);
  fprintf( Param .output_param, "Ds is %s\n",              Param .Ds);
#ifndef NI
  fprintf (Param.output_param, "a is %lg\n",               Param .a);
  fprintf (Param.output_param, "k is %lg\n",               Param .k);
#endif
  fprintf (Param.output_param, "StoreInterPos is %lg\n",   Param .StoreInterPos);
  fprintf (Param.output_param, "StoreInterDisp is %lg\n",  Param .StoreInterDisp);
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
  fprintf(Param.output_pos, "\n%.4f",State._time);
  for (n=0;n<Param.N;n++){
    fprintf(Param.output_pos,"\t%lg",Particles[n].dx);
  }
  fprintf(Param.output_pos,"\n%.4f",State._time);
  for (n=0;n<Param.N;n++){
    fprintf(Param.output_pos,"\t%lg",Particles[n].dy);
  }
  fprintf(Param.output_pos,"\n");
  fflush(Param.output_pos);
}



// Store the identities and integrated displacements of all particles, with the time
void Store_Disp(param Param, particle* Particles, state State){
  long n;
  fprintf(Param.output_disp,"%.4f",State._time);
  for (n=0;n<Param.N;n++) fprintf(Param.output_disp,"\t%lg",Particles[n].idx);
  fprintf(Param.output_disp,"\n%.4f",State._time);
  for (n=0;n<Param.N;n++) fprintf(Param.output_disp,"\t%lg",Particles[n].idy);
  fprintf(Param.output_disp,"\n");
  fflush(Param.output_disp);
}


void Update_Profile(param Param, particle* Particles, data* Data){
  long i; // particle index
  int j,k;// index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    Data[0].prof[j][k] += 1;
  }
}

void Store_Profile(param Param,particle* Particles, state State, data* Data){
  int j,k; // profile index
  
  for (j=0; j<Param.Nbinx; j++){
    fprintf(Param.output_prof, "%.4f",State._time);
    for (k=0; k<Param.Nbiny; k++){
      fprintf(Param.output_prof, "\t%ld", Data[0].prof[j][k]);
      Data[0].prof[j][k]=0;
    }
    fprintf(Param.output_prof, "\n");
  }
  fflush(Param.output_prof);
}


void Update_J(param Param, particle* Particles, data* Data){
  long i; // particle index
  int j,k;// index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
#ifndef NI
    Data[0].Jprof[j][k][0] += Particles[i].ifx;
    Data[0].Jprof[j][k][1] += Particles[i].ify;
#endif
    Data[0].Jprof[j][k][2] += Particles[i].vx*Param.dt;
    Data[0].Jprof[j][k][3] += Particles[i].vy*Param.dt;
  }
}

void Store_J(param Param,particle* Particles, state State, data* Data){
  int n,j,k; // profile index
  int n0 = 0; // which J piece to start with
#ifdef NI 
  n0 = 2;
#endif
  
  for (n=n0; n<4; n++){
    for (j=0; j<Param.Nbinx; j++){
      fprintf(Param.output_Jprof, "%.4f",State._time);
      for (k=0; k<Param.Nbiny; k++){
        fprintf(Param.output_Jprof, "\t%lg", Data[0].Jprof[j][k][n]);
        Data[0].Jprof[j][k][n]=0.0;
      }
      fprintf(Param.output_Jprof, "\n");
    }
  }
  fflush(Param.output_Jprof);
}



void Update_Data(param Param, particle* Particles, state* State, data* Data){
  if (State[0]._time > State[0].NextUpdateProf-EPS){
    Update_Profile(Param, Particles, Data);
    Update_J(Param, Particles, Data);
    State[0].NextUpdateProf += Param.UpdateInterProf;
  }
}

void Store_Data(param Param, particle* Particles, state* State, data* Data){
  if (State[0]._time > State[0].NextStorePos - EPS){
    Store_Pos(Param, Particles, State[0]);
    State[0].NextStorePos += Param.StoreInterPos;
  }

  if (State[0]._time > State[0].NextStoreDisp - EPS){
    Store_Disp(Param, Particles, State[0]);
    State[0].NextStoreDisp += Param.StoreInterDisp;
  }

  if (State[0]._time > State[0].NextStoreProf-EPS){
    Store_Profile(Param, Particles, State[0], Data);
    Store_J(Param, Particles, State[0], Data);
    State[0].NextStoreProf += Param.StoreInterProf;
  }
}


/* Print percentage of completed simulation */
void Update_Progress(param Param, state* State){
  printf("\rIn progress [%d %%]", (int)(State[0].prev_percentage + 1.01));
  fflush(stdout);
  State[0].prev_percentage += 1.0;
}