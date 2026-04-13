
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
#if defined(TOFX_LIN) || defined(TOFX_CUB)
  strcat(command_base, "Nstep_T ");                  argctarget ++;
  strcat(command_base, "T_xs ");                     argctarget ++;
  strcat(command_base, "Ts ");                       argctarget ++;
#elif defined(TOFX_TRIG)
  strcat(command_base, "Nmode_T ");                  argctarget ++;
  strcat(command_base, "T0 ");                       argctarget ++;
  strcat(command_base, "T_a10,T_a01,T_a20,... ");    argctarget ++;
  strcat(command_base, "T_bx10,T_bx01,T_bx20,... "); argctarget ++;
  strcat(command_base, "T_by10,T_by01,T_by20,... "); argctarget ++;
#else
  strcat(command_base, "T0 ");                       argctarget ++;
#endif
#if defined(FOFX_LIN) || defined(FOFX_CUB)
  strcat(command_base, "Nstep_F ");                  argctarget ++;
  strcat(command_base, "F_xs ");                     argctarget ++;
  strcat(command_base, "Fs ");                       argctarget ++;
#endif
#ifndef NI
  strcat (command_base, "a ");                       argctarget ++;
  strcat (command_base, "k ");                       argctarget ++;
#endif
  strcat (command_base, "StoreInterPos ");           argctarget ++;
#ifdef STORE_DISP
  strcat(command_base, "StoreInterDisp ");           argctarget ++;
#endif
#if defined(STORE_ANYPROF)
  strcat(command_base, "Nbinx ");                    argctarget ++;
  strcat(command_base, "Nbiny ");                    argctarget ++;
  strcat(command_base, "NstepProf ");                argctarget ++;
  strcat(command_base, "StoreInterProf ");           argctarget ++;
  strcat(command_base, "UpdateInterProf ");          argctarget ++;
#endif
#ifdef TRACK
  strcat(command_base, "Ntrack ");                   argctarget ++;
#endif
#ifdef GIVEN_IC
  strcat(command_base, "IC_file ");                  argctarget ++;
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
  Param[0].output_param = fopen (name, "w");

  sprintf (name, "%s-pos", argv[i]);
  Param[0].output_pos = fopen (name, "w");

#ifdef STORE_DISP
  sprintf (name, "%s-disp", argv[i]);
  Param[0].output_disp = fopen (name, "w");
#endif

#ifdef STORE_PROF
  sprintf (name, "%s-prof", argv[i]);
  Param[0].output_prof = fopen (name, "w");
#endif

#ifdef STORE_JPROF
  sprintf (name, "%s-Jprof", argv[i]);
  Param[0].output_Jprof = fopen (name, "w");
#endif

#ifdef STORE_FPROF
  sprintf (name, "%s-Fprof", argv[i]);
  Param[0].output_Fprof = fopen (name, "w");
#endif

#ifdef STORE_TPROF
  sprintf (name, "%s-Tprof", argv[i]);
  Param[0].output_Tprof = fopen (name, "w");
#endif
  
#ifdef STORE_SIGMAIKPROF
  sprintf(name,"%s-sigmaIKprof",argv[i]);
  Param[0].output_sigmaIKprof=fopen(name,"w");
#endif
  
#ifdef STORE_EPR
  sprintf(name,"%s-EPR",argv[i]);
  Param[0].output_EPR=fopen(name,"w");
#endif

#ifdef TRACK
  sprintf (name, "%s-traj", argv[i]);
  Param[0].output_traj = fopen (name, "w");
#endif

  printf ("created files\n");

  i++;

  // Read in basic parameters
  Param[0].Lx                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Ly                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].N                  = (long) strtod (argv[i], NULL)      ; i++ ;
  Param[0].dt                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].tf                 = strtod (argv[i], NULL)             ; i++ ;
  Param[0].seed               = (long long) strtod (argv[i], NULL) ; i++ ;
#if defined(TOFX_LIN) || defined(TOFX_CUB)
  Param[0].TParam.Nstep       = (int) strtod(argv[i], NULL)        ; i++ ;
  Param[0].Txs                = argv[i]                            ; i++ ;
  Param[0].Ts                 = argv[i]                            ; i++ ;
#elif defined(TOFX_TRIG)
  Param[0].TParam.Nmode       = (int) strtod(argv[i], NULL)        ; i++ ;
  Param[0].TParam.f0          = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Tas                = argv[i]                            ; i++ ;
  Param[0].Tbxs               = argv[i]                            ; i++ ;
  Param[0].Tbys               = argv[i]                            ; i++ ;
#else
  Param[0].T0                 = strtod (argv[i], NULL)             ; i++ ;
#endif
#if defined(FOFX_LIN) || defined(FOFX_CUB)
  Param[0].FParam.Nstep       = (int) strtod(argv[i], NULL)        ; i++ ;
  Param[0].Fxs                = argv[i]                            ; i++ ;
  Param[0].Fs                 = argv[i]                            ; i++ ;
#endif
#ifndef NI
  Param[0].a                  = strtod (argv[i], NULL)             ; i++ ;
  Param[0].k                  = strtod (argv[i], NULL)             ; i++ ;
#endif
  Param[0].StoreInterPos      = strtod (argv[i], NULL)             ; i++ ;
#ifdef STORE_DISP
  Param[0].StoreInterDisp     = strtod (argv[i], NULL)             ; i++ ;
#endif
#if defined(STORE_ANYPROF)
  Param[0].Nbinx              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].Nbiny              = strtod (argv[i], NULL)             ; i++ ;
  Param[0].NstepProf          = (int)strtod (argv[i], NULL)        ; i++ ;
  Param[0].StoreInterProf     = strtod (argv[i], NULL)             ; i++ ;
  Param[0].UpdateInterProf    = strtod (argv[i], NULL)             ; i++ ;
#endif
#ifdef TRACK
  Param[0].Ntrack             = strtod (argv[i], NULL)             ; i++ ;
#endif
#ifdef GIVEN_IC
  Param[0].IC_file            = argv[i]                            ; i++ ;
#endif
  

  // define parameters that are functions of the inputs and/or known

  Particles[0]             = (particle *)malloc (sizeof (particle) * Param[0].N);
  init_genrand64 (Param[0].seed);
  State[0]._time           = 0;
  State[0].prev_percentage = 0;
  State[0].NextStorePos    = Param[0].StoreInterPos;
  Param[0].UpdateInterProg = Param[0].tf/100.0; // update every 1%
  State[0].NextUpdateProg  = Param[0].UpdateInterProg;
#ifdef STORE_DISP
  State[0].NextStoreDisp   = Param[0].StoreInterDisp;
#endif

#if defined(STORE_ANYPROF)
  State[0].NextStoreProf   = Param[0].StoreInterProf;
  State[0].NextUpdateProf  = State[0].NextStoreProf - (Param[0].UpdateInterProf) * (Param[0].NstepProf - 1);
  assert((Param[0].NstepProf-1)*Param[0].UpdateInterProf <= Param[0].StoreInterProf);
#endif

#ifdef TOFX 
  Param[0].sqrt2dt              = pow(2*Param[0].dt, 0.5);
#endif
#if defined(TOFX_LIN) || defined(TOFX_CUB)
  Param[0].TParam.L             = Param[0].Lx;
  Param[0].TParam.xs            = (double*) calloc(Param[0].TParam.Nstep, sizeof(double));
  Param[0].TParam.fs            = (double*) calloc(Param[0].TParam.Nstep, sizeof(double));
  split_vec(Param[0].TParam.Nstep, &Param[0].TParam.xs, Param[0].Txs, ',');
  split_vec(Param[0].TParam.Nstep, &Param[0].TParam.fs, Param[0].Ts, ',');
#endif
#ifdef TOFX_LIN
  compute_coefs_lin(&Param[0].TParam);
#elif defined(TOFX_CUB)
  compute_coefs_cub(&Param[0].TParam);
#elif defined(TOFX_TRIG)
  Param[0].TParam.Lx            = Param[0].Lx;
  Param[0].TParam.Ly            = Param[0].Ly;
  Param[0].TParam.amps          = (double*) calloc(Param[0].TParam.Nmode, sizeof(double));
  Param[0].TParam.phases_x      = (double*) calloc(Param[0].TParam.Nmode, sizeof(double));
  Param[0].TParam.phases_y      = (double*) calloc(Param[0].TParam.Nmode, sizeof(double));
  split_vec(Param[0].TParam.Nmode, &Param[0].TParam.amps, Param[0].Tas, ',');
  split_vec(Param[0].TParam.Nmode, &Param[0].TParam.phases_x, Param[0].Tbxs, ',');
  split_vec(Param[0].TParam.Nmode, &Param[0].TParam.phases_y, Param[0].Tbys, ',');
  compute_coefs_trig(&Param[0].TParam);
#else
  Param[0].sqrt2T0dt = pow(2 * Param[0].T0 * Param[0].dt, 0.5);
#endif




#if defined(FOFX_LIN) || defined(FOFX_CUB)
  Param[0].FParam.L             = Param[0].Lx;
  Param[0].FParam.xs            = (double*) calloc(Param[0].FParam.Nstep, sizeof(double));
  Param[0].FParam.fs            = (double*) calloc(Param[0].FParam.Nstep, sizeof(double));
  split_vec(Param[0].FParam.Nstep, &Param[0].FParam.xs, Param[0].Fxs, ',');
  split_vec(Param[0].FParam.Nstep, &Param[0].FParam.fs, Param[0].Fs, ',');
  for (i=0; i<Param[0].FParam.Nstep; i++){ // scale by dt
    Param[0].FParam.fs[i] = Param[0].FParam.fs[i] * Param[0].dt;
  }
#endif
#ifdef FOFX_LIN
  compute_coefs_lin(&Param[0].FParam);
#elif defined(FOFX_CUB)
  compute_coefs_cub(&Param[0].FParam);
#endif


#ifndef NI
  
#ifdef FORCE_WCA
  Param[0].rmax            = Param[0].a * pow(2, 1/6.);
#else
  Param[0].rmax            = Param[0].a;
#endif
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

#ifdef FORCE_WCA
  get_force_WCA_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, 
                       &Param[0].c1, &Param[0].c2);
#elif defined(FORCE_QUARTIC)
  get_force_quartic_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, 
                           &Param[0].c1, &Param[0].c2, &Param[0].c3);
#elif defined(FORCE_CUBIC)
  get_force_cubic_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, 
                         &Param[0].c1, &Param[0].c2, &Param[0].c3);
#elif defined(FORCE_HARMONIC)
  get_force_harmonic_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, 
                            &Param[0].c1, &Param[0].c2);
#elif defined(FORCE_IHARMONIC)
  get_force_iharmonic_params(Param[0].k, Param[0].a, 1.0, Param[0].dt, &Param[0].c1);
#elif defined(FORCE_LINEAR)
  get_force_linear_params(Param[0].k, 1.0, Param[0].dt, &Param[0].c1);
#endif
#endif


#if defined(STORE_PROF) || defined(STORE_JPROF)
  Param[0].binwidthx = (double)Param[0].Lx / ((double)Param[0].Nbinx);
  Param[0].binwidthy = (double)Param[0].Ly / ((double)Param[0].Nbiny);
#endif

#if defined(STORE_PROF)
  Data[0].prof = calloc (Param[0].Nbinx, sizeof (long **));
  for (i = 0; i < Param[0].Nbinx; i++) {
    Data[0].prof[i] = calloc (Param[0].Nbiny, sizeof (long *));
    for (j = 0; j < Param[0].Nbiny; j++) {
      Data[0].prof[i][j] = 0.0;
    }
  }
#endif
#if defined(STORE_JPROF)
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
#endif
#if defined(STORE_FPROF)
  Data[0].Fprof = calloc(Param[0].Nbinx, sizeof(double**));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].Fprof[i] = calloc(Param[0].Nbiny, sizeof(double*));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].Fprof[i][j] = calloc(2, sizeof(double));
      Data[0].Fprof[i][j][0] = 0.0; // fx
      Data[0].Fprof[i][j][1] = 0.0; // fy
    }
  }
#endif
#if defined(STORE_TPROF)
  Data[0].Tprof = calloc(Param[0].Nbinx, sizeof(double*));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].Tprof[i] = calloc(Param[0].Nbiny, sizeof(double));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].Tprof[i][j] = 0.0;
    }
  }
#endif

#if defined(STORE_SIGMAIKPROF)
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
#endif

#if defined(STORE_EPR)
  Data[0].EPR = calloc(Param[0].Nbinx, sizeof(double*));
  for (i=0; i<Param[0].Nbinx; i++){
    Data[0].EPR[i] = calloc(Param[0].Nbiny, sizeof(double));
    for (j=0; j<Param[0].Nbiny; j++){
      Data[0].EPR[i][j] = 0.0;
    }
  }
#endif


// initialize particle locations
#ifdef RANDOM_IC
  for (i = 0; i < Param[0].N; i++) {

      // random initialize positions uniformly over interval
      Particles[0][i].x = Param[0].Lx * genrand64_real2 ();
      Particles[0][i].y = Param[0].Ly * genrand64_real2 ();

#if !defined(TOFX_LIN) & !defined(TOFX_CUB) & !defined(TOFX_TRIG)
      Particles[0][i].T = Param[0].T0; // temperature
#endif
      Particles[0][i].dx = 0.0;
      Particles[0][i].dy = 0.0;
#ifdef STORE_DISP
      Particles[0][i].idx = 0.0;
      Particles[0][i].idy = 0.0;
#endif
#if defined(STORE_JPROF)
      Particles[0][i].dx_prev = 0.0;
      Particles[0][i].dy_prev = 0.0;
#endif
#ifdef FOFX
      Particles[0][i].Fx = 0.0;
      Particles[0][i].Fy = 0.0;
#endif

#ifndef NI
      // Add particle to its box and link to its neighbors
      Particles[0][i].bi = (int)(floor (Particles[0][i].x / Param[0].rmax) + EPS);
      Particles[0][i].bj = (int)(floor (Particles[0][i].y / Param[0].rmax) + EPS);
      AddinBox (i, Particles[0][i].bi, Particles[0][i].bj, Boxes[0], Neighbours[0]);

      Particles[0][i].ifx = 0.0;
      Particles[0][i].ify = 0.0;
#ifdef STORE_SIGMAIKPROF
      Particles[0][i].ifxs = (double*) calloc(Param[0].N, sizeof(double));
      Particles[0][i].ifys = (double*) calloc(Param[0].N, sizeof(double));
      for (j=0; j<Param[0].N; j++) {
        Particles[0][i].ifxs[j] = 0.0;
        Particles[0][i].ifys[j] = 0.0;
      }
#endif
#endif
    }
#endif
}

#if defined(TOFX_LIN)
double Tofx(fparam_lin TParam, double x, double y){
  return fofx_lin(TParam, x);
}
#elif defined(TOFX_CUB)
double Tofx(fparam_cub TParam, double x, double y){
  return fofx_cub(TParam, x);
}
#elif defined(TOFX_TRIG)
double Tofx(fparam_trig TParam, double x, double y){
  return fofx_trig(TParam, x, y);
}
#endif

#if defined(FOFX_LIN)
double Fofx_dt(fparam_lin FParam, double x, double y){
  return fofx_lin(FParam, x);
}
#elif defined(FOFX_CUB)
double Fofx_dt(fparam_cub FParam, double x, double y){
  return fofx_cub(FParam, x);
}
#endif


void Update_Particles(param Param,particle* Particles, state State
#ifndef NI
                      , long* Neighbours, long** Boxes
#endif
                      ){
  long n;
  int newbi,newbj;
  double T,sqrt2Tdt;

  // reset displacements and forces to 0
  for (n=0; n<Param.N; n++){
#if defined(STORE_JPROF)
    Particles[n].dx_prev = Particles[n].dx;
    Particles[n].dy_prev = Particles[n].dy;
#endif
    Particles[n].dx = 0.0;
    Particles[n].dy = 0.0;
#ifndef NI
    Particles[n].ifx = 0.0;
    Particles[n].ify = 0.0;
#ifdef STORE_SIGMAIKPROF
    memset(Particles[n].ifxs, 0.0, Param.N);
    memset(Particles[n].ifys, 0.0, Param.N);
#endif
#endif
  }
  
  // incorporate interactions
#ifndef NI
  Loop_Force(Param, Particles, Neighbours, Boxes);
#endif

  for (n=0 ; n<Param.N ; n++ ){

    // Compute the displacement of each particle
#if defined(TOFX)
    Particles[n].T = Tofx(Param.TParam, Particles[n].x, Particles[n].y);
    sqrt2Tdt = pow(Particles[n].T, 0.5) * Param.sqrt2dt;
#else
    sqrt2Tdt = Param.sqrt2T0dt;
#endif
    Particles[n].dx += sqrt2Tdt * gasdev();
    Particles[n].dy += sqrt2Tdt * gasdev();

#ifdef FOFX
    // incorporate external force
    Particles[n].Fx = Fofx_dt(Param.FParam, Particles[n].x, Particles[n].y);
    Particles[n].dx += Particles[n].Fx;
#endif

#ifndef NI
    // incorporate interactions
    Particles[n].dx += Particles[n].ifx;
    Particles[n].dy += Particles[n].ify;
#endif
  }
  
  for ( n=0 ; n<Param.N ; n++ ){
    // Once all displacement are computed, move the particles
    Particles[n].x += Particles[n].dx;
    Particles[n].y += Particles[n].dy;

    // Take care of periodic boundary conditions
#ifndef CBC_X
    while(Particles[n].x>Param.Lx) Particles[n].x-=Param.Lx;
    while(Particles[n].x<0) Particles[n].x+=Param.Lx;
#endif
#ifndef CBC_Y
    while(Particles[n].y>Param.Ly) Particles[n].y-=Param.Ly;
    while(Particles[n].y<0) Particles[n].y+=Param.Ly;
#endif

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

#ifdef STORE_DISP
    // Update the cumulated displacements
    Particles[n].idx += Particles[n].dx;
    Particles[n].idy += Particles[n].dy;
#endif
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
#if defined(TOFX_LIN) || defined(TOFX_CUB)
  fprintf( Param .output_param, "Nstep_T is %d\n",         Param .TParam .Nstep);
  fprintf( Param .output_param, "Txs is %s\n",             Param .Txs);
  fprintf( Param .output_param, "Ts is %s\n",              Param .Ts);
#elif defined(TOFX_TRIG)
  fprintf (Param.output_param, "Nmode_T is %d\n",          Param .TParam.Nmode);
  fprintf (Param.output_param, "T0 is %lg\n",              Param .TParam.f0);
  fprintf (Param.output_param, "T_as is %s\n",             Param .Tas);
  fprintf (Param.output_param, "T_bxs is %s\n",            Param .Tbxs);
  fprintf (Param.output_param, "T_bys is %s\n",            Param .Tbys);
#else
  fprintf (Param.output_param, "T0 is %lg\n",              Param .T0);
#endif
#if defined(FOFX_LIN) || defined(FOFX_CUB)
  fprintf( Param .output_param, "Nstep_F is %d\n",         Param .FParam .Nstep);
  fprintf( Param .output_param, "Fxs is %s\n",             Param .Fxs);
  fprintf( Param .output_param, "Fs is %s\n",              Param .Fs);
#endif
#ifndef NI
  fprintf (Param.output_param, "a is %lg\n",               Param .a);
  fprintf (Param.output_param, "k is %lg\n",               Param .k);
#endif
  fprintf (Param.output_param, "StoreInterPos is %lg\n",   Param .StoreInterPos);
#ifdef STORE_DISP
  fprintf (Param.output_param, "StoreInterDisp is %lg\n",  Param .StoreInterDisp);
#endif
#if defined(STORE_ANYPROF)
  fprintf (Param.output_param, "Nbinx is %d\n",            Param .Nbinx);
  fprintf (Param.output_param, "Nbiny is %d\n",            Param .Nbiny);
  fprintf (Param.output_param, "NstepProf is %d\n",        Param .NstepProf);
  fprintf (Param.output_param, "StoreInterProf is %lg\n",  Param .StoreInterProf);
  fprintf (Param.output_param, "UpdateInterProf is %lg\n", Param .UpdateInterProf);
#endif
#ifdef TRACK
  fprintf (Param.output_param, "Ntrack is %lg\n",          Param .Ntrack);
#endif
#ifdef GIVEN_IC
  fprintf (Param.output_param, "IC_file is %lg\n",         Param .IC_file);
#endif
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


#ifdef TRACK
  // Store the complete trajectories (position and polarization) of first Ntrack
  // particles
void Store_Trajectories (param Param, particle * Particles, state State){
  long n; // particle index

  for (n = 0; n < Param.Ntrack; n++) {
    fprintf (Param.output_traj, "%lg\t%ld\t%lg\t%lg\n", State._time, n, Particles[n].x,
             Particles[n].y);
    }
}
#endif


#ifdef STORE_DISP
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
#endif

#if defined(STORE_PROF)
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
#endif

#ifdef STORE_JPROF
void Update_J(param Param, particle* Particles, data* Data){
  long i; // particle index
  int j,k;// index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  double Jxdt,Jydt;
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    Jxdt = Particles[i].dx_prev - (Particles[i].dx-Particles[i].dx_prev)/Param.dt;
    Jydt = Particles[i].dy_prev - (Particles[i].dy-Particles[i].dy_prev)/Param.dt;
#ifndef NI
    Data[0].Jprof[j][k][0] += Particles[i].ifx;
    Data[0].Jprof[j][k][1] += Particles[i].ify;
    Data[0].Jprof[j][k][2] += Jxdt - Particles[i].ifx;
    Data[0].Jprof[j][k][3] += Jydt - Particles[i].ify;
#else
    Data[0].Jprof[j][k][2] += Jxdt;
    Data[0].Jprof[j][k][3] += Jydt;
#endif
  }
}

void Store_J(param Param,particle* Particles, state State, data* Data){
  int n,j,k; // profile index
  
  for (n=0; n<4; n++){
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
#endif


#ifdef STORE_FPROF
void Update_F(param Param, particle* Particles, data* Data){
  long i;
  int j,k;
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    Data[0].Fprof[j][k][0] += Particles[i].Fx;
    Data[0].Fprof[j][k][1] += Particles[i].Fy;
  }
}
void Store_F(param Param,particle* Particles, state State, data* Data){
  int n,j,k; // profile index
  
  for (n=0; n<2; n++){
    for (j=0; j<Param.Nbinx; j++){
      fprintf(Param.output_Fprof, "%.4f",State._time);
      for (k=0; k<Param.Nbiny; k++){
        fprintf(Param.output_Fprof, "\t%lg", Data[0].Fprof[j][k][n]);
        Data[0].Fprof[j][k][n]=0.0;
      }
      fprintf(Param.output_Fprof, "\n");
    }
  }
  fflush(Param.output_Fprof);
}
#endif

#ifdef STORE_TPROF
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
#endif


#ifdef STORE_SIGMAIKPROF
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
#endif



#ifdef STORE_EPR 

// calculate EPR -- requires Stratonovich average, so moving particles back by dx/2
void Update_EPR(param Param, particle* Particles, data* Data
#ifndef NI
                      , long* Neighbours, long** Boxes
#endif
  ){
  long i;  // particle index
  int j,k; // profile index
  int newbi,newbj; // box membership
  double fx, fy; // total x and y force

  // move to position at middle of timestep
  for (i=0; i<Param.N; i++){
    Particles[i].x -= Particles[i].dx * 0.5;
    Particles[i].y -= Particles[i].dy * 0.5;

    Particles[i].ifx = 0.0;
    Particles[i].ify = 0.0;

    // Take care of periodic boundary conditions
#ifndef CBC_X
    while(Particles[i].x>Param.Lx) Particles[i].x-=Param.Lx;
    while(Particles[i].x<0) Particles[i].x+=Param.Lx;
#endif
#ifndef CBC_Y
    while(Particles[i].y>Param.Ly) Particles[i].y-=Param.Ly;
    while(Particles[i].y<0) Particles[i].y+=Param.Ly;
#endif


#ifndef NI 
    // Update box membership
    newbi = (int) (floor(Particles[i].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[i].y/Param.rmax) + EPS);

    if(Particles[i].bi!=newbi || Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Boxes,Neighbours);
      AddinBox(i,newbi,newbj,Boxes,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
#endif
  }

#ifndef NI 
  // calculate interaction force at middle of timestep
  Loop_Force(Param, Particles, Neighbours, Boxes);
#endif

  // calculate EPR
  for (i=0; i<Param.N; i++){
    fx = Particles[i].ifx;
    fy = Particles[i].ify;

#ifdef FOFX 
    // add external x force at intermediate location
    fx += Fofx_dt(Param.FParam, Particles[i].x, Particles[i].y);
#endif

#ifdef TOFX 
    // calculate temperature at intermediate location
    Particles[i].T = Tofx(Param.TParam, Particles[i].x, Particles[i].y);
#endif

    j = (int) (floor(Particles[i].x/Param.binwidthx) + EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy) + EPS);
    Data[0].EPR[j][k] += (Particles[i].dx * fx + Particles[i].dy * fy) / (Particles[i].T * Param.dt);
  }

  // move back to original position 
  for (i=0; i<Param.N; i++){
    Particles[i].x += Particles[i].dx * 0.5;
    Particles[i].y += Particles[i].dy * 0.5;

    // Take care of periodic boundary conditions
#ifndef CBC_X
    while(Particles[i].x>Param.Lx) Particles[i].x-=Param.Lx;
    while(Particles[i].x<0) Particles[i].x+=Param.Lx;
#endif
#ifndef CBC_Y
    while(Particles[i].y>Param.Ly) Particles[i].y-=Param.Ly;
    while(Particles[i].y<0) Particles[i].y+=Param.Ly;
#endif

#ifndef NI 
    // Update box membership
    newbi = (int) (floor(Particles[i].x/Param.rmax) + EPS);
    newbj = (int) (floor(Particles[i].y/Param.rmax) + EPS);

    if(Particles[i].bi!=newbi || Particles[i].bj!=newbj){
      RemovefromBox(i,Particles[i].bi,Particles[i].bj,Boxes,Neighbours);
      AddinBox(i,newbi,newbj,Boxes,Neighbours);
      Particles[i].bi=newbi;
      Particles[i].bj=newbj;
    }
#endif
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
#endif



void Update_Data(param Param, particle* Particles, state* State, data* Data
#if !defined(NI) & defined(STORE_EPR)
                      , long* Neighbours, long** Boxes
#endif
  ){
#if defined(STORE_ANYPROF)
  if (State[0]._time > State[0].NextUpdateProf-EPS){
#ifdef STORE_PROF
    Update_Profile(Param, Particles, Data);
#endif
#ifdef STORE_JPROF
    Update_J(Param, Particles, Data);
#endif
#ifdef STORE_FPROF
    Update_F(Param, Particles, Data);
#endif
#ifdef STORE_TPROF
    Update_T(Param, Particles, Data);
#endif
#ifdef STORE_SIGMAIKPROF
    Update_SigmaIK(Param, Particles, Data);
#endif
#ifdef STORE_EPR 
    Update_EPR(Param, Particles, Data
#if !defined(NI)
              , Neighbours, Boxes
#endif
      );
#endif
    State[0].NextUpdateProf += Param.UpdateInterProf;
  }
#endif
}

void Store_Data(param Param, particle* Particles, state* State, data* Data){
  if (State[0]._time > State[0].NextStorePos - EPS){
    Store_Pos(Param, Particles, State[0]);
    State[0].NextStorePos += Param.StoreInterPos;
  }

#ifdef STORE_DISP 
  if (State[0]._time > State[0].NextStoreDisp - EPS){
    Store_Disp(Param, Particles, State[0]);
    State[0].NextStoreDisp += Param.StoreInterDisp;
  }
#endif

#if defined(STORE_ANYPROF)
  if (State[0]._time > State[0].NextStoreProf-EPS){
#ifdef STORE_PROF 
    Store_Profile(Param, Particles, State[0], Data);
#endif
#ifdef STORE_JPROF
    Store_J(Param, Particles, State[0], Data);
#endif
#ifdef STORE_FPROF
    Store_F(Param, Particles, State[0], Data);
#endif
#ifdef STORE_TPROF
    Store_T(Param, Particles, State[0], Data);
#endif
#ifdef STORE_SIGMAIKPROF
    Store_SigmaIK(Param, Particles, State[0], Data);
#endif
#ifdef STORE_EPR
    Store_EPR(Param, Particles, State[0], Data);
#endif
    State[0].NextStoreProf += Param.StoreInterProf;
  }
#endif

#ifdef TRACK 
  Store_Trajectories(Param, Particles, State[0]);
#endif
}


/* Print percentage of completed simulation */
void Update_Progress(param Param, state* State){
  printf("\rIn progress [%d %%]", (int)(State[0].prev_percentage + 1.01));
  fflush(stdout);
  State[0].prev_percentage += 1.0;
}