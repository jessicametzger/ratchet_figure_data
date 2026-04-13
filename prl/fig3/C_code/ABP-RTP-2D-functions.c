
void Initialize_Parameters(int argc, char* argv[], param* Param, state* State, data* Data, particle** Particles, 
  long*** Boxes, long** Neighbours
  ){
  
  long i,j;
  double tmp;
  int argctarget=0;           // Number of parameters that should be used
  char command_base[1000]=""; // string that contains the desired format of the command line
  char name[200];             // string in which the file names are written

  /*
    Format the command line and count the arguments
  */
  argctarget = 0;
  strcat(command_base, "usage: ");
  strcat(command_base, argv[0]);                  argctarget ++;
  strcat(command_base, " ");
  strcat(command_base, "file ");                  argctarget ++;
  strcat(command_base, "Lx ");                    argctarget ++;
  strcat(command_base, "Ly ");                    argctarget ++;
  strcat(command_base, "N ");                     argctarget ++;
  strcat(command_base, "dt ");                    argctarget ++;
  strcat(command_base, "tf ");                    argctarget ++;
  strcat(command_base, "tau ");                   argctarget ++;
  strcat(command_base, "Nstep_v ");               argctarget ++;
  strcat(command_base, "v_x1,v_x2,... ");         argctarget ++;
  strcat(command_base, "v1,v2,... ");             argctarget ++;
  strcat(command_base, "a ");                     argctarget ++;
  strcat(command_base, "k ");                     argctarget ++;
  strcat(command_base, "seed ");                  argctarget ++;
  strcat(command_base, "StoreInterPos ");         argctarget ++;
  strcat(command_base, "Nbinx ");                 argctarget ++;
  strcat(command_base, "Nbiny ");                 argctarget ++;
  strcat(command_base, "NstepProf ");             argctarget ++;
  strcat(command_base, "StoreInterProf ");        argctarget ++;
  strcat(command_base, "UpdateInterProf ");       argctarget ++;
  strcat(command_base, "StoreInterDisp ");        argctarget ++;

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
  Param[0].output_param=fopen(name,"w");

  sprintf(name,"%s-pos",argv[i]);
  Param[0].output_pos=fopen(name,"w");
  
  sprintf(name,"%s-prof",argv[i]);
  Param[0].output_prof=fopen(name,"w");

  sprintf(name,"%s-FAprof",argv[i]);
  Param[0].output_FAprof=fopen(name,"w");

  sprintf(name,"%s-Fintprof",argv[i]);
  Param[0].output_Fintprof=fopen(name,"w");

  sprintf(name,"%s-sigmaAprof",argv[i]);
  Param[0].output_sigmaAprof=fopen(name,"w");

  sprintf(name,"%s-sigmaIKprof",argv[i]);
  Param[0].output_sigmaIKprof=fopen(name,"w");

  sprintf(name,"%s-disp",argv[i]);
  Param[0].output_disp=fopen(name,"w");
  
  printf("created files\n");
  
  i++;
  
  // Read in basic parameters
  Param[0].Lx                 = strtod(argv[i], NULL)        ; i ++;
  Param[0].Ly                 = strtod(argv[i], NULL)        ; i ++;
  Param[0].N                  = (long) strtod(argv[i], NULL) ; i ++;
  Param[0].dt                 = strtod(argv[i], NULL)        ; i ++;
  Param[0].tf                 = strtod(argv[i], NULL)        ; i ++;
  Param[0].tau                = strtod(argv[i], NULL)        ; i ++;
  Param[0].vParam.Nstep       = (int) strtod(argv[i], NULL)  ; i ++;
  Param[0].vxs                = argv[i]                      ; i ++;
  Param[0].vs                 = argv[i]                      ; i ++;
  Param[0].seed               = strtod(argv[i], NULL)        ; i ++;
  Param[0].StoreInterPos      = strtod(argv[i], NULL)        ; i ++;
  Param[0].Nbinx              = (int) strtod(argv[i], NULL)  ; i ++;
  Param[0].Nbiny              = (int) strtod(argv[i], NULL)  ; i ++;
  Param[0].NstepProf          = (int) strtod(argv[i], NULL)  ; i ++;
  Param[0].StoreInterProf     = strtod(argv[i], NULL)        ; i ++;
  Param[0].UpdateInterProf    = strtod(argv[i], NULL)        ; i ++;
  Param[0].StoreInterDisp     = strtod(argv[i], NULL)        ; i ++;

  
  /*
    Define parameters that are functions of inputs and/or known
  */

  init_genrand64(Param[0].seed);
  State[0]._time                 = 0.0;
  State[0].prev_percentage       = 0.0;
  State[0].NextStorePos          = Param[0].StoreInterPos;
  Param[0].UpdateInterProg       = Param[0].tf/100.0; // update every 1%
  State[0].NextUpdateProg        = Param[0].UpdateInterProg;
  Particles[0]                   = (particle*) malloc(sizeof(particle)*Param[0].N);
  Param[0].two_pi                = 2*M_PI;

  State[0].NextStoreDisp         = Param[0].StoreInterDisp;

  Param[0].binwidthx       = (double) Param[0].Lx/((double) Param[0].Nbinx);
  Param[0].binwidthy       = (double) Param[0].Ly/((double) Param[0].Nbiny);
  State[0].NextStoreProf      = Param[0].StoreInterProf;
  State[0].NextUpdateProf     = State[0].NextStoreProf - Param[0].UpdateInterProf * (Param[0].NstepProf-1);
  
  // Must have enough time to make enough measurements between profile storages
  assert(Param[0].NstepProf*Param[0].UpdateInterProf<=Param[0].StoreInterProf);

  Param[0].sqrt2Drdt             = sqrt(2*Param[0].dt/Param[0].tau);


  Param[0].vParam.L             = Param[0].Lx;
  Param[0].vParam.xs            = (double*) calloc(Param[0].vParam.Nstep, sizeof(double));
  Param[0].vParam.fs            = (double*) calloc(Param[0].vParam.Nstep, sizeof(double));
  split_vec(Param[0].vParam.Nstep, &Param[0].vParam.xs, Param[0].vxs, ',');
  split_vec(Param[0].vParam.Nstep, &Param[0].vParam.fs, Param[0].vs, ',');
  for (i=0; i<Param[0].vParam.Nstep; i++){
    Param[0].vParam.fs[i]          = Param[0].vParam.fs[i]*Param[0].dt; // scale by dt
  }
  Param[0].vParam.cs            = (double*) calloc(Param[0].vParam.Nstep, sizeof(double));
  Param[0].vParam.ds            = (double*) calloc(Param[0].vParam.Nstep, sizeof(double));
  for (i=0; i<Param[0].vParam.Nstep-1; i++){
    tmp                         = Param[0].vParam.xs[i+1] - Param[0].vParam.xs[i];
    Param[0].vParam.cs[i]       = -3 * (Param[0].vParam.fs[i] - Param[0].vParam.fs[i+1]) / pow(tmp,2);
    Param[0].vParam.ds[i]       = 2 * (Param[0].vParam.fs[i] - Param[0].vParam.fs[i+1]) / pow(tmp,3);
  }
  tmp                           = Param[0].vParam.xs[0] + Param[0].Lx - Param[0].vParam.xs[i];
  Param[0].vParam.ds[i]         = 2 * (Param[0].vParam.fs[i] - Param[0].vParam.fs[0]) / pow(tmp,3);
  Param[0].vParam.cs[i]         = -3 * (Param[0].vParam.fs[i] - Param[0].vParam.fs[0]) / pow(tmp,2);


  
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
  
  
  Data[0].prof              = (long**) malloc(Param[0].Nbinx*sizeof(long*));
  for (i=0;i<Param[0].Nbinx;i++){
    Data[0].prof[i]         = (long*) calloc(Param[0].Nbiny,sizeof(long));
    for (j=0;j<Param[0].Nbiny;j++){
      Data[0].prof[i][j]    = 0;
    }
  }

  Data[0].FAprof            = (double***) malloc(Param[0].Nbinx*sizeof(double**));
  for (i=0;i<Param[0].Nbinx;i++){
    Data[0].FAprof[i]       = (double**) calloc(Param[0].Nbiny,sizeof(double*));
    for (j=0;j<Param[0].Nbiny;j++){
      Data[0].FAprof[i][j]  = (double*) calloc(2, sizeof(double));
      Data[0].FAprof[i][j][0] = 0.0;  // ux*(xdot*dvdx + ydot*dvdy)
      Data[0].FAprof[i][j][1] = 0.0;  // uy*(xdot*dvdx + ydot*dvdy)
    }
  }

  Data[0].Fintprof            = (double***) malloc(Param[0].Nbinx*sizeof(double**));
  for (i=0;i<Param[0].Nbinx;i++){
    Data[0].Fintprof[i]       = (double**) calloc(Param[0].Nbiny,sizeof(double*));
    for (j=0;j<Param[0].Nbiny;j++){
      Data[0].Fintprof[i][j]  = (double*) calloc(2, sizeof(double));
      Data[0].Fintprof[i][j][0] = 0.0;  // ifx
      Data[0].Fintprof[i][j][1] = 0.0;  // ify
    }
  }

  Data[0].sigmaAprof            = (double***) malloc(Param[0].Nbinx*sizeof(double**));
  for (i=0;i<Param[0].Nbinx;i++){
    Data[0].sigmaAprof[i]       = (double**) calloc(Param[0].Nbiny,sizeof(double*));
    for (j=0;j<Param[0].Nbiny;j++){
      Data[0].sigmaAprof[i][j]  = (double*) calloc(4, sizeof(double));
      Data[0].sigmaAprof[i][j][0] = 0.0;  // ux*xdot*v*tau
      Data[0].sigmaAprof[i][j][1] = 0.0;  // ux*ydot*v*tau
      Data[0].sigmaAprof[i][j][2] = 0.0;  // uy*xdot*v*tau
      Data[0].sigmaAprof[i][j][3] = 0.0;  // uy*ydot*v*tau
    }
  }
  
  Data[0].sigmaIKprof              = (double***) malloc(Param[0].NxBox*sizeof(double**));
  for (i=0;i<Param[0].NxBox;i++){
    Data[0].sigmaIKprof[i]         = (double**) calloc(Param[0].NyBox,sizeof(double*));
    for (j=0;j<Param[0].NyBox;j++){
      Data[0].sigmaIKprof[i][j]    = (double*) calloc(3, sizeof(double));
      Data[0].sigmaIKprof[i][j][0] = 0.0; // sigmaIK_xx
      Data[0].sigmaIKprof[i][j][1] = 0.0; // sigmaIK_xy
      Data[0].sigmaIKprof[i][j][2] = 0.0; // sigmaIK_yy
    }
  }


// initialize particle locations and polarizations  
  for (i=0;i<Param[0].N;i++){
  
    // random initialize positions uniformly over interval
    Particles[0][i].x=Param[0].Lx*genrand64_real2();
    Particles[0][i].y=Param[0].Ly*genrand64_real2();
    
    // randomly pick direction
    Particles[0][i].th = 2*M_PI*genrand64_real2();

    Particles[0][i].v_dt = 0.0;
    Particles[0][i].ux = cos(Particles[0][i].th);
    Particles[0][i].uy = sin(Particles[0][i].th);
    Particles[0][i].dx = 0.0;
    Particles[0][i].dy = 0.0;
    Particles[0][i].ifx = 0.0;
    Particles[0][i].ify = 0.0;
    Particles[0][i].ifxs = (double*) calloc(Param[0].N,sizeof(double));
    Particles[0][i].ifys = (double*) calloc(Param[0].N,sizeof(double));
    for (j=0; j<Param[0].N; j++){
      Particles[0][i].ifxs[j] = 0.0;
      Particles[0][i].ifys[j] = 0.0;
    }

    Particles[0][i].idx = 0.0;
    Particles[0][i].idy = 0.0;

    // Add particle to its box and link to its neighbors
    Particles[0][i].bi= (int) (floor(Particles[0][i].x/Param[0].rmax) + EPS);
    Particles[0][i].bj= (int) (floor(Particles[0][i].y/Param[0].rmax) + EPS);
    AddinBox(i,Particles[0][i].bi,Particles[0][i].bj,Boxes[0],Neighbours[0]);
  }
}



double vofx_dt(fparam_cub vParam, double x, double y){
  return fofx_cub(vParam, x);
}
double dvdx_dt(fparam_cub vParam, double x, double y){
  return dfdx_cub(vParam, x);
}
double dvdy_dt(fparam_cub vParam, double x, double y){
  return 0.0;
}



void Update_Particles(param Param,particle* Particles, state State, long* Neighbours, long** Boxes){
  long n;
  int newbi,newbj;
  double v_dt;

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
    Particles[n].dx = 0.0;
    Particles[n].dy = 0.0;
  
    /* Change angle of self-propulsion */
    Particles[n].th += Param.sqrt2Drdt * gasdev();
    
    if(Particles[n].th>Param.two_pi) Particles[n].th -= Param.two_pi;
    if(Particles[n].th<0)            Particles[n].th += Param.two_pi;

    Particles[n].ux = cos(Particles[n].th);
    Particles[n].uy = sin(Particles[n].th);

    // Compute the displacement of each particle as x += v cos(theta) dt, y += v sin(theta) dt
    Particles[n].v_dt = vofx_dt(Param.vParam, Particles[n].x, Particles[n].y);
    Particles[n].dx += Particles[n].v_dt * Particles[n].ux;
    Particles[n].dy += Particles[n].v_dt * Particles[n].uy;

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

    // Update the cumulated displacements
    Particles[n].idx += Particles[n].dx;
    Particles[n].idy += Particles[n].dy;
  }
}

void Store_Parameters(int argc, char* argv[], param Param){
  long i;
  
  /* Store parameters */
  for(i=0;i<argc;i++){
    fprintf(Param.output_param,"%s ",argv[i]);
  }
  fprintf( Param .output_param,"\n");

  fprintf( Param .output_param,"Lx is %lg\n",                 Param .Lx);
  fprintf( Param .output_param,"Ly is %lg\n",                 Param .Ly);
  fprintf( Param .output_param,"N is %ld\n",                  Param .N);
  fprintf( Param .output_param,"dt is %lg\n",                 Param .dt);
  fprintf( Param .output_param,"tf is %lg\n",                 Param .tf);
  fprintf( Param .output_param,"tau is %lg\n",                Param .tau);
  fprintf( Param .output_param, "Nstep_v is %d\n",            Param .vParam .Nstep);
  fprintf( Param .output_param, "vx1,vx2,... is %s\n",        Param .vxs);
  fprintf( Param .output_param, "v1,v2,... is %s\n",          Param .vs);
  fprintf( Param .output_param,"a is %lg\n",                  Param .a);
  fprintf( Param .output_param,"k is %lg\n",                  Param .k);
  fprintf( Param .output_param,"seed is %lld\n",              Param .seed);
  fprintf( Param .output_param,"StoreInterPos is %lg\n",      Param .StoreInterPos);
  fprintf( Param .output_param,"Nbinx is %d\n",               Param .Nbinx);
  fprintf( Param .output_param,"Nbiny is %d\n",               Param .Nbiny);
  fprintf( Param .output_param,"NstepProf is %d\n",           Param .NstepProf);
  fprintf( Param .output_param,"StoreInterProf is %lg\n",     Param .StoreInterProf);
  fprintf( Param .output_param,"UpdateInterProf is %lg\n",    Param .UpdateInterProf);
  fprintf( Param .output_param,"StoreInterDisp is %lg\n",     Param .StoreInterDisp);
  fflush(  Param .output_param);
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
    Data[0].prof[j][k]     += 1;
  }
}
void Store_Profile(param Param,state State, data* Data){
  int n,j,k; // profile index
  
  /* density profile */
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


void Update_Fint(param Param, particle* Particles, data* Data){
  long i;  // particle index
  int j,k; // index to add to

  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    Data[0].Fintprof[j][k][0] += Particles[i].ifx;
    Data[0].Fintprof[j][k][1] += Particles[i].ify;
  }
}
void Store_Fint(param Param, state State, data* Data){
  int n,j,k; // profile index

  for (n=0; n<2; n++){
    for (j=0; j<Param.Nbinx; j++){
      fprintf(Param.output_Fintprof,"%.4f",State._time);
      for (k=0; k<Param.Nbiny; k++){
        fprintf(Param.output_Fintprof, "\t%lg", Data[0].Fintprof[j][k][n]);
        Data[0].Fintprof[j][k][n] = 0.0;
      }
      fprintf(Param.output_Fintprof,"\n");
    }
  }
  fflush(Param.output_Fintprof);
}


void Update_FA(param Param, particle* Particles, data* Data){
  long i;  // particle index
  int j,k; // index to add to
  double dvdx,dvdy; // slope of v

  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    dvdx = dvdx_dt(Param.vParam, Particles[i].x, Particles[i].y)/Param.dt;
    dvdy = dvdy_dt(Param.vParam, Particles[i].x, Particles[i].y)/Param.dt;
    Data[0].FAprof[j][k][0] += Particles[i].ux*Particles[i].dx*dvdx + Particles[i].ux*Particles[i].dy*dvdy;
    Data[0].FAprof[j][k][1] += Particles[i].uy*Particles[i].dx*dvdx + Particles[i].uy*Particles[i].dy*dvdy;
  }
}
void Store_FA(param Param,state State, data* Data){
  int n,j,k; // profile index
  
  for (n=0; n<2; n++){
    for (j=0; j<Param.Nbinx; j++){
      fprintf(Param.output_FAprof, "%.4f",State._time);
      for (k=0; k<Param.Nbiny; k++){
        fprintf(Param.output_FAprof, "\t%lg", Data[0].FAprof[j][k][n]);
        Data[0].FAprof[j][k][n]=0.0;
      }
      fprintf(Param.output_FAprof, "\n");
    }
  }
  fflush(Param.output_FAprof);
}


void Update_SigmaA(param Param, particle* Particles, data* Data){
  long i; // particle index
  int j,k;// index to add to. From 0 to Param.Nbinx-1 and 0 to Param.Nbiny-1
  double v_dt_tau;
  
  for (i=0; i<Param.N; i++){
    j = (int) (floor(Particles[i].x/Param.binwidthx)+EPS);
    k = (int) (floor(Particles[i].y/Param.binwidthy)+EPS);
    v_dt_tau = Param.tau*Particles[i].v_dt;
    Data[0].sigmaAprof[j][k][0] += Particles[i].ux*Particles[i].dx*v_dt_tau;
    Data[0].sigmaAprof[j][k][1] += Particles[i].ux*Particles[i].dy*v_dt_tau;
    Data[0].sigmaAprof[j][k][2] += Particles[i].uy*Particles[i].dx*v_dt_tau;
    Data[0].sigmaAprof[j][k][3] += Particles[i].uy*Particles[i].dy*v_dt_tau;
  }
}

void Store_SigmaA(param Param,state State, data* Data){
  int n,j,k; // profile index
  
  for (n=0; n<4; n++){
    for (j=0; j<Param.Nbinx; j++){
      fprintf(Param.output_sigmaAprof, "%.4f",State._time);
      for (k=0; k<Param.Nbiny; k++){
        fprintf(Param.output_sigmaAprof, "\t%lg", Data[0].sigmaAprof[j][k][n]);
        Data[0].sigmaAprof[j][k][n]=0.0;
      }
      fprintf(Param.output_sigmaAprof, "\n");
    }
  }
  fflush(Param.output_sigmaAprof);
}



void Store_SigmaIK(param Param,state State, data* Data) {
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


void Update_Data(param Param, particle* Particles, state* State, data* Data){

  if (State[0]._time > State[0].NextUpdateProf-EPS){
    Update_Profile(Param, Particles, Data);
    Update_FA(Param, Particles, Data);
    Update_Fint(Param, Particles, Data);
    Update_SigmaA(Param, Particles, Data);
    Update_SigmaIK(Param, Particles, Data);
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
    Store_Profile(Param, State[0], Data);
    Store_FA(Param, State[0], Data);
    Store_Fint(Param, State[0], Data);
    Store_SigmaA(Param, State[0], Data);
    Store_SigmaIK(Param, State[0], Data);
    State[0].NextStoreProf += Param.StoreInterProf;
  }
}


/* Print percentage of completed simulation */
void Update_Progress(param Param, state* State){
  printf("\rIn progress [%d %%]", (int)(State[0].prev_percentage + 1.01));
  fflush(stdout);
  State[0].prev_percentage += 1.0;
}