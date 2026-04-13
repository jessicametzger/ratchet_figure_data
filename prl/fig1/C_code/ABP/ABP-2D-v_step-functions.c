
void Initialize_Parameters(int argc, char* argv[], param* Param, state* State, data* Data, particle** Particles
#ifndef NI
  , long*** Boxes, long** Neighbours
#endif
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
#ifndef NI
  strcat(command_base, "a ");                     argctarget ++;
  strcat(command_base, "k ");                     argctarget ++;
#endif
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
#ifndef NI
  Param[0].a                  = strtod(argv[i], NULL)        ; i ++;
  Param[0].k                  = strtod(argv[i], NULL)        ; i ++;
#endif
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
  Param[0].vParam.dfdxs         = (double*) calloc(Param[0].vParam.Nstep, sizeof(double));
  for (i=0; i<Param[0].vParam.Nstep-1; i++){
    Param[0].vParam.dfdxs[i]    = (Param[0].vParam.fs[i+1] - Param[0].vParam.fs[i]) / (Param[0].vParam.xs[i+1] - Param[0].vParam.xs[i]);
  }
  Param[0].vParam.dfdxs[i]      = (Param[0].vParam.fs[0] - Param[0].vParam.fs[i]) / (Param[0].vParam.xs[0]+Param[0].Lx - Param[0].vParam.xs[i]);


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
  
  Data[0].prof              = (long**) malloc(Param[0].Nbinx*sizeof(long*));
  for (i=0;i<Param[0].Nbinx;i++){
    Data[0].prof[i]         = (long*) calloc(Param[0].Nbiny,sizeof(long));
    for (j=0;j<Param[0].Nbiny;j++){
      Data[0].prof[i][j]    = 0;
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
#ifndef NI
    Particles[0][i].ifx = 0.0;
    Particles[0][i].ify = 0.0;
#endif
    Particles[0][i].idx = 0.0;
    Particles[0][i].idy = 0.0;
    
#ifndef NI
    // Add particle to its box and link to its neighbors
    Particles[0][i].bi= (int) (floor(Particles[0][i].x/Param[0].rmax) + EPS);
    Particles[0][i].bj= (int) (floor(Particles[0][i].y/Param[0].rmax) + EPS);
    AddinBox(i,Particles[0][i].bi,Particles[0][i].bj,Boxes[0],Neighbours[0]);
#endif
  }
}


double vofx_dt(fparam_lin vParam, double x, double y){
  return fofx_lin(vParam, x);
}


void Update_Particles(param Param,particle* Particles, state State
#ifndef NI
                      , long* Neighbours, long** Boxes
#endif
                      ){
  long n;
  int newbi,newbj;
  double v_dt;

  // reset displacements and forces to 0
  for (n=0; n<Param.N; n++){
    Particles[n].dx = 0.0;
    Particles[n].dy = 0.0;
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
#ifndef NI
  fprintf( Param .output_param,"a is %lg\n",                  Param .a);
  fprintf( Param .output_param,"k is %lg\n",                  Param .k);
#endif
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


void Update_Data(param Param, particle* Particles, state* State, data* Data){
  if (State[0]._time > State[0].NextUpdateProf-EPS){
    Update_Profile(Param, Particles, Data);
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
    State[0].NextStoreProf += Param.StoreInterProf;
  }
}


/* Print percentage of completed simulation */
void Update_Progress(param Param, state* State){
  printf("\rIn progress [%d %%]", (int)(State[0].prev_percentage + 1.01));
  fflush(stdout);
  State[0].prev_percentage += 1.0;
}