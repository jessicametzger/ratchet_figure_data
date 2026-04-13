/*
	Pairwise forces functions.
*/


/*
  WCA: F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6) + epsilon
  
  so the force magnitude is
  
  F = 4 epsilon (12 sigma^12 / r^13 - 6 sigma^6 / r^7)
  
  We use the LJ potential with a cut-off at r=sigma*2^1/6
*/
#ifdef FORCE_WCA
// c1=4*12*dt*amp*sigma^12, c2=-24*amp*sigma^6
double get_force_WCA_params(double k, double a, double mu, double dt, double* c1, double* c2){
  c1[0] = 4*12*dt*k*mu*pow(a,12);
  c2[0] = -24*dt*k*mu*pow(a,6);
}
double Force_WCA(double dist, double dist2, double c1, double c2){
  double dist6=dist2*dist2*dist2;
  double dist7=dist6*dist;
  double dist13=dist7*dist6;
  return c1/dist13 + c2/dist7;
}
#endif


/*
  QUARTIC:  F = -grad V where

  V(r) = k*(-3(r/a)^4 + 4(r/a)^3 (1+r0/a) - 6(r/a)^2 (r0/a) + 2r0/a - 1) / (2r0/a - 1)

  so that F(0)=F(r0)=F(a)=0, V(0) = 1, and V(a) = 0.
  If r0=a/2, throws error. If r0<a/2, then it's a local maximum. 
  If r0>a, then there's no local minimum

  see: https://www.desmos.com/calculator/2eptglk46v
*/
#ifdef FORCE_QUARTIC 
double get_force_quartic_params(double k, double a, double mu, double dt, double r0,
  double* c1, double* c2, double* c3){
  c3[0] = k * dt * mu * 12 / pow(a,4) / (2*r0/a - 1);
  c2[0] = k * dt * mu * (-12) * (1+r0/a) / pow(a,3) / (2*r0/a - 1);
  c1[0] = k * dt * mu * 12 * (r0/a) / pow(a,2) / (2*r0/a - 1);
}
double Force_quartic(double dist, double dist2, double c1, double c2, double c3){
  double dist3 = dist2*dist;
	return c3*dist3 + c2*dist2 + c1*dist;
}
#endif


/*
  CUBIC:  F = -grad V where

  V(r) = k*(-2(r/a)^3 + 3(r/a)^2 (1+r0/a) - 6(r/a)(r0/a) + 3r0/a - 1) / (3r0/a - 1)

  so that F(r0)=F(a)=0, V(0) = 1, and V(a) = 0.

  If r0=a/3, throws error. If r0<a/3, then it's a local maximum. 
  If r0>a, then there's no local minimum

  see: https://www.desmos.com/calculator/2eptglk46v
*/
#ifdef FORCE_CUBIC
double get_force_cubic_params(double k, double a, double mu, double dt, double r0,
  double* c1, double* c2, double* c3){
  c3[0] = k * dt * mu * 6 / pow(a,3) / (3*r0/a - 1);
  c2[0] = k * dt * mu * (-6) * (1+r0/a) / pow(a,2) / (3*r0/a - 1);
  c1[0] = k * dt * mu * 6 * (r0/a) / a / (3*r0/a - 1);
}
double Force_cubic(double dist, double dist2, double c1, double c2, double c3){
  double dist3 = dist2*dist;
  return c3*dist3 + c2*dist2 + c1*dist;
}
#endif


/*
  HARMONIC: We use harmonic spheres F = -\grad V where
  
  V(r)= (amp sigma / 2) (1 - r/sigma)^2

  and force magnitude

  F(r_i-r_j) = amp (1 - |r_i-r_j|/sigma)
  
  We use a cut-off at r=sigma 
*/
#ifdef FORCE_HARMONIC
double get_force_harmonic_params(double k, double a, double mu, double dt, double* c1, double* c2){
  c1[0] = k * dt * mu;
  c2[0] = -k * dt * mu / a;
}
double Force_harmonic(double dist, double c1, double c2){
	return c1 + c2*dist;
}
#endif



/*
  INV_HARMONIC: We use V(r) = (amp sigma / 2) (1 - (r/sigma)^2)
  
  and force magnitude
  
  F(r_i-r_j) = amp |r_i-r_j| / sigma
  
  We use a cut-off at r=sigma
*/
#ifdef FORCE_IHARMONIC
double get_force_iharmonic_params(double k, double a, double mu, double dt, double* c1){
  c1[0] = k * dt * mu / a;
}
double Force_iharmonic(double dist, double c1){
	return c1*dist;
}
#endif


/*
  "tent" potential
  
  V(r) = amp sigma (1-|r|/sigma)
  
  F(r) = amp 
  
  We use a cut-off at r=sigma
*/
#ifdef FORCE_LINEAR
double get_force_linear_params(double k, double mu, double dt, double* c1){
  c1[0] = k * dt * mu;
}
double Force_linear(double c1){
	return c1;
}
#endif


#ifndef ONED

// given (corrected) distances dx,dy and indices j,k, update forces (stored in "Particles").
// ASSUME THEY ARE CLOSE ENOUGH TO INTERACT!
void Pair_force(double dx, double dy, double dist2, long j, long k, particle* Particles, param Param){
  double Force,Forcex,Forcey;

  double dist=pow(dist2,0.5);

#ifdef FORCE_WCA
  Force = Force_WCA(dist,dist2,Param.c1,Param.c2);
#elif defined(FORCE_QUARTIC)
  Force = Force_quartic(dist,dist2,Param.c1,Param.c2,Param.c3);
#elif defined(FORCE_CUBIC)
  Force = Force_cubic(dist,dist2,Param.c1,Param.c2,Param.c3);
#elif defined(FORCE_IHARMONIC)
  Force = Force_iharmonic(dist,Param.c1);
#elif defined(FORCE_HARMONIC)
  Force = Force_harmonic(dist,Param.c1,Param.c2);
#elif defined(FORCE_LINEAR)
  Force = Force_linear(Param.c1);
#endif
  Forcex = Force*dx/dist;
  Forcey = Force*dy/dist;
#ifdef STORE_SIGMAIKPROF
  Particles[j].ifxs[k] = Forcex;
  Particles[j].ifys[k] = Forcey;
  Particles[k].ifxs[j] = -Forcex;
  Particles[k].ifys[j] = -Forcey;
#endif
  Particles[j].ifx    += Forcex;
  Particles[j].ify    += Forcey;
  Particles[k].ifx    -= Forcex;
  Particles[k].ify    -= Forcey;
}


void Compute_force_neighbors(param Param, particle* Particles, long j,int bi,int bj,
                             long* Neighbours, long** Boxes){
  int m;                      //neighbour box id
  int nbi,nbj;                //box name
  double dx_box,dy_box;       // particle offset
  long k;                     //particle index
  double Force,Forcex,Forcey; // Force exerted between two particles
  double dx,dy,dist2;         // displacement between particles
  
  // NeighbouringBoxes[bi][bj][0] is box (bi,bj) itself
  for (m=0;m<5;m++){
    nbi=Param.NeighbouringBoxes[bi][bj][m].i;           // x index of the mth neighbor
    nbj=Param.NeighbouringBoxes[bi][bj][m].j;           // y index of the mth neighbor
    dx_box=Param.NeighbouringBoxes[bi][bj][m].epsilonx; // x offset to be added to the particles in
                                                        // box (bi,bj) to take into account the periodic
                                                        // boundary conditions
    dy_box=Param.NeighbouringBoxes[bi][bj][m].epsilony; // y offset to be added to the particles in
                                                        // box (bi,bj) to take into account the periodic
                                                        // boundary conditions

    //Loop through all particles k in box 
    if (m==0) k = Neighbours[2*j+1]; // Start with next in same box
    else      k = Boxes[nbi][nbj];     //Start with the first in neighbouring box
  
    //As long as there are particles in the box, iterate
    while(k!=-1){

      dx = Particles[j].x - Particles[k].x - dx_box; // positive if x_j>x_k
      dy = Particles[j].y - Particles[k].y - dy_box; // positive if y_j>y_k

      //If the particles are closer than Param.rmax, they interact
      if ((m!=0) & ((dx>Param.rmax) | (dy>Param.rmax) | (dx<-Param.rmax) | dy<-Param.rmax)) {
        k=Neighbours[2*k+1]; // k is now the next particle in the list
        continue;
      }

      dist2=pow(dx,2) + pow(dy,2);
      if (dist2>Param.rmax2){
        k=Neighbours[2*k+1]; // k is now the next particle in the list
        continue;
      }

      // calculate the force
      Pair_force(dx,dy,dist2,j,k,Particles,Param);
      k=Neighbours[2*k+1]; // k is now the next particle in the list
    }
  }
}


void Loop_Force(param Param, particle* Particles, long* Neighbours, long** Boxes){
  int bi,bj; // Box index
  long j; // Particle index

  for (bi=0; bi<Param.NxBox; bi++){
    for (bj=0; bj<Param.NyBox; bj++){

      //Loop through all the particles 'j' in the box.
      j=Boxes[bi][bj]; //Start with j being the 1st particle

      //As long as j is not -1, compute its interactions with all the other particles
      while(j!=-1){
      
        // Loop through the boxes
        Compute_force_neighbors(Param, Particles, j, bi, bj, Neighbours, Boxes);
        
        //We have included all the interactions between j and other particles. Now iterate over j
        j=Neighbours[2*j+1];
      }
      // We are done with all particles in Boxes[bi][bj], iterate over the box
    }
  }
}


#ifdef STORE_SIGMAIKPROF
void Update_SigmaIK(param Param, particle* Particles, data* Data){
  long j,k;                                                     // Particle ids
  int bi,bj,nbi,nbj,nbi1,nbj1,m;                                // box coordinates
  double x,y,xn,yn;                                             // coordinates (self & neighbour)
  double sigma_xx, sigma_xy, sigma_yy;                          // stress tensor components
  double frac_0, frac_2, frac_4;                                // frac_k = fraction of separation vector in box k
  double sigma_xx_part, sigma_xy_part, sigma_yy_part;           // stress tensor components (piece within box)
  double sigma_xx_part1, sigma_xy_part1, sigma_yy_part1;        // stress tensor components (piece within neighbor box)
  double box_left, box_bottom, box_right, box_top, intersect_x; // coordinates of the box sides
  double dx,dy,dist2,dist,dx_box,dy_box;                        // displacements and PBC s

  for (j=0; j<Param.N; j++){

#if !defined(CENTER_X) & !defined(CENTER_Y)
    x = Particles[j].x - Particles[j].dx;
    y = Particles[j].y - Particles[j].dy;
#else
    x = Particles[j].xshift - Particles[j].dx;
    y = Particles[j].yshift - Particles[j].dy;
#endif
#ifndef CBC_X
    while(x>Param.Lx) x-=Param.Lx;
    while(x<0) x+=Param.Lx;
#endif
#ifndef CBC_Y
    while(y>Param.Ly) y-=Param.Ly;
    while(y<0) y+=Param.Ly;
#endif
    bi = (int) (floor(x/Param.rmax) + EPS);
    bj = (int) (floor(y/Param.rmax) + EPS);
    
    box_left = Param.rmax * floor(x/Param.rmax);
    box_bottom = Param.rmax * floor(y/Param.rmax);
    box_right = box_left + Param.rmax;
    box_top = box_bottom + Param.rmax;

    for (k=0; k<Param.N; k++){
      if (k==j) continue;
#if !defined(CENTER_X) & !defined(CENTER_Y)
      xn = Particles[k].x - Particles[k].dx;
      yn = Particles[k].y - Particles[k].dy;
#else
      xn = Particles[k].xshift - Particles[k].dx;
      yn = Particles[k].yshift - Particles[k].dy;
#endif
#ifndef CBC_X
      while(xn>Param.Lx) xn-=Param.Lx;
      while(xn<0) xn+=Param.Lx;
#endif
#ifndef CBC_Y
      while(yn>Param.Ly) yn-=Param.Ly;
      while(yn<0) yn+=Param.Ly;
#endif
      nbi = (int) (floor(xn/Param.rmax) + EPS);
      nbj = (int) (floor(yn/Param.rmax) + EPS);

      for (m=0; m<5; m++){
        if ((nbi==Param.NeighbouringBoxes[bi][bj][m].i) & (nbj==Param.NeighbouringBoxes[bi][bj][m].j))
          break;
      }
      if (m==5) continue; // if k isn't a neighbor of j
      if ((m==0) & (k>j)) continue; // remove duplicates

      dx_box = Param.NeighbouringBoxes[bi][bj][m].epsilonx;
      dy_box = Param.NeighbouringBoxes[bi][bj][m].epsilony;

      dx=x-dx_box-xn;// positive if x_j>x_k
      dy=y-dy_box-yn;// positive if y_j>y_k

      dist2=pow(dx,2) + pow(dy,2);
      if (dist2>Param.rmax2) continue; // if k doesn't interact with j

      sigma_xx = dx*Particles[j].ifxs[k];
      sigma_xy = dx*Particles[j].ifys[k];
      sigma_yy = dy*Particles[j].ifys[k];


      if (m==0){ // same box
        Data[0].sigmaIKprof[bi][bj][0] += sigma_xx;
        Data[0].sigmaIKprof[bi][bj][1] += sigma_xy;
        Data[0].sigmaIKprof[bi][bj][2] += sigma_yy;

      } else if (m==1){ //above
        // dy will be negative
        frac_0 = (box_top-y)/(-dy);

        // add to same box
        sigma_xx_part = sigma_xx*frac_0;
        sigma_xy_part = sigma_xy*frac_0;
        sigma_yy_part = sigma_yy*frac_0;
        
        Data[0].sigmaIKprof[bi][bj][0] += sigma_xx_part;
        Data[0].sigmaIKprof[bi][bj][1] += sigma_xy_part;
        Data[0].sigmaIKprof[bi][bj][2] += sigma_yy_part;

        // add to neighbor box
        Data[0].sigmaIKprof[nbi][nbj][0] += sigma_xx-sigma_xx_part;
        Data[0].sigmaIKprof[nbi][nbj][1] += sigma_xy-sigma_xy_part;
        Data[0].sigmaIKprof[nbi][nbj][2] += sigma_yy-sigma_yy_part;

      } else if (m==2){ //above and to the right
        // dx and dy will be negative
        intersect_x = (box_top - y)*(-dx)/(-dy);
        if (intersect_x > box_right){ // it goes thru box 3
          frac_0 = (box_right-x)/(-dx);
          frac_2 = (yn+dy_box-box_top)/(-dy);
          nbi1 = Param.NeighbouringBoxes[bi][bj][3].i; // box 3
          nbj1 = Param.NeighbouringBoxes[bi][bj][3].j; // box 3
        } else { // goes thru box 1
          frac_0 = (box_top-y)/(-dy);
          frac_2 = (xn+dx_box-box_right)/(-dx);
          nbi1 = Param.NeighbouringBoxes[bi][bj][1].i; // box 1
          nbj1 = Param.NeighbouringBoxes[bi][bj][1].j; // box 1
        }

        // add components (same box)
        sigma_xx_part = sigma_xx*frac_0;
        sigma_xy_part = sigma_xy*frac_0;
        sigma_yy_part = sigma_yy*frac_0;
        
        Data[0].sigmaIKprof[bi][bj][0] += sigma_xx_part;
        Data[0].sigmaIKprof[bi][bj][1] += sigma_xy_part;
        Data[0].sigmaIKprof[bi][bj][2] += sigma_yy_part;

        // add to neighbor box
        sigma_xx_part1 = sigma_xx*frac_2;
        sigma_xy_part1 = sigma_xy*frac_2;
        sigma_yy_part1 = sigma_yy*frac_2;
        Data[0].sigmaIKprof[nbi][nbj][0] += sigma_xx_part1;
        Data[0].sigmaIKprof[nbi][nbj][1] += sigma_xy_part1;
        Data[0].sigmaIKprof[nbi][nbj][2] += sigma_yy_part1;

        // add to overlapping box
        Data[0].sigmaIKprof[nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
        Data[0].sigmaIKprof[nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
        Data[0].sigmaIKprof[nbi1][nbj1][2] += sigma_yy-sigma_yy_part-sigma_yy_part1;

      } else if (m==3){ //to the right
        // dx will be negative
        frac_0 = (box_right-x)/(-dx);

        // add to same box
        sigma_xx_part = sigma_xx*frac_0;
        sigma_xy_part = sigma_xy*frac_0;
        sigma_yy_part = sigma_yy*frac_0;
        
        Data[0].sigmaIKprof[bi][bj][0] += sigma_xx_part;
        Data[0].sigmaIKprof[bi][bj][1] += sigma_xy_part;
        Data[0].sigmaIKprof[bi][bj][2] += sigma_yy_part;

        // add to neighbor box
        Data[0].sigmaIKprof[nbi][nbj][0] += sigma_xx-sigma_xx_part;
        Data[0].sigmaIKprof[nbi][nbj][1] += sigma_xy-sigma_xy_part;
        Data[0].sigmaIKprof[nbi][nbj][2] += sigma_yy-sigma_yy_part;

      } else if (m==4){ //below and to the right
        // dx will be negative and dy will be positive
        intersect_x = (y - box_bottom)*(-dx)/dy;

        if (intersect_x < box_right) { // goes thru box "5" (box below)
          frac_0 = (y - box_bottom)/dy; // how much goes to self
          frac_4 = (xn + dx_box - box_right)/(-dx); // how much goes to neighbor (box 4)

          // will overlap the box below.
          nbi1 = bi;
          nbj1 = (bj-1+Param.NyBox)%Param.NyBox;

        } else { // goes thru box 3
          frac_0 = (box_right - x)/(-dx); // how much goes to self
          frac_4 = (box_bottom - (yn+dy_box))/dy; // how much goes to neighbor (box 4)
        }

        // add components (same box)
        sigma_xx_part = sigma_xx*frac_0;
        sigma_xy_part = sigma_xy*frac_0;
        sigma_yy_part = sigma_yy*frac_0;
        
        Data[0].sigmaIKprof[bi][bj][0] += sigma_xx_part;
        Data[0].sigmaIKprof[bi][bj][1] += sigma_xy_part;
        Data[0].sigmaIKprof[bi][bj][2] += sigma_yy_part;

        // add to neighbor box
        sigma_xx_part1 = sigma_xx*frac_4;
        sigma_xy_part1 = sigma_xy*frac_4;
        sigma_yy_part1 = sigma_yy*frac_4;
        Data[0].sigmaIKprof[nbi][nbj][0] += sigma_xx_part1;
        Data[0].sigmaIKprof[nbi][nbj][1] += sigma_xy_part1;
        Data[0].sigmaIKprof[nbi][nbj][2] += sigma_yy_part1;

        // add to overlapping box
        Data[0].sigmaIKprof[nbi1][nbj1][0] += sigma_xx-sigma_xx_part-sigma_xx_part1;
        Data[0].sigmaIKprof[nbi1][nbj1][1] += sigma_xy-sigma_xy_part-sigma_xy_part1;
        Data[0].sigmaIKprof[nbi1][nbj1][2] += sigma_yy-sigma_yy_part-sigma_yy_part1;
      }
    }
  }
}
#endif

#elif defined(ONED)

// given (corrected) distance dx and indices j,k, update forces (stored in "Particles").
// ASSUME THEY ARE CLOSE ENOUGH TO INTERACT!
void Pair_force(double dx, double dist, long j, long k, particle* Particles, param Param){
  double Force,Forcex,Forcey;

#if defined(FORCE_CUBIC) || defined(FORCE_QUARTIC) || defined(FORCE_WCA)
  double dist2=dist*dist;
#endif

#ifdef FORCE_WCA
  Force = Force_WCA(dist,dist2,Param.c1,Param.c2);
#elif defined(FORCE_QUARTIC)
  Force = Force_quartic(dist,dist2,Param.c1,Param.c2,Param.c3);
#elif defined(FORCE_CUBIC)
  Force = Force_cubic(dist,dist2,Param.c1,Param.c2,Param.c3);
#elif defined(FORCE_IHARMONIC)
  Force = Force_iharmonic(dist,Param.c1);
#elif defined(FORCE_HARMONIC)
  Force = Force_harmonic(dist,Param.c1,Param.c2);
#elif defined(FORCE_LINEAR)
  Force = Force_linear(Param.c1);
#endif
  Forcex = (dx>0)?Force:(-Force);
#ifdef STORE_SIGMAIKPROF
  Particles[j].ifxs[k] = Forcex;
  Particles[k].ifxs[j] = -Forcex;
#endif
  Particles[j].ifx    += Forcex;
  Particles[k].ifx    -= Forcex;
}


// spatial hashing
#ifndef NO_SH

void Compute_force_neighbors(param Param, particle* Particles, long j,int bi,
                             long* Neighbours, long* Boxes){
  int m;                      //neighbour box id
  int nbi;                    //box name
  double dx_box;              // particle offset
  long k;                     //particle index
  double dx,dist;             // displacement between particles
  
  // NeighbouringBoxes[bi][bj][0] is box (bi,bj) itself
  for (m=0;m<2;m++){
    nbi=Param.NeighbouringBoxes[bi][m].i;           // x index of the mth neighbor
    dx_box=Param.NeighbouringBoxes[bi][m].epsilonx; // x offset to be added to the particles in

    //Loop through all particles k in box 
    if (m==0) k = Neighbours[2*j+1]; // Start with next in same box
    else      k = Boxes[nbi];        //Start with the first in neighbouring box
  
    //As long as there are particles in the box, iterate
    while(k!=-1){

      dx = Particles[j].x - Particles[k].x - dx_box; // positive if x_j>x_k

      //If the particles are closer than Param.rmax, they interact
      dist = (dx>=0)?dx:(-dx);
      if (dist>=Param.rmax){
        k=Neighbours[2*k+1]; // k is now the next particle in the list
        continue;
      }

      // calculate the force
      Pair_force(dx,dist,j,k,Particles,Param);
      k=Neighbours[2*k+1]; // k is now the next particle in the list
    }
  }
}


void Loop_Force(param Param, particle* Particles, long* Neighbours, long* Boxes){
  int bi; // Box index
  long j; // Particle index

  for (bi=0; bi<Param.NxBox; bi++){

    //Loop through all the particles 'j' in the box.
    j=Boxes[bi]; //Start with j being the 1st particle

    //As long as j is not -1, compute its interactions with all the other particles
    while(j!=-1){
    
      // Loop through the boxes
      Compute_force_neighbors(Param, Particles, j, bi, Neighbours, Boxes);
      
      //We have included all the interactions between j and other particles. Now iterate over j
      j=Neighbours[2*j+1];
    }
    // We are done with all particles in Boxes[bi][bj], iterate over the box
  }
}

#else // no spatial hashing
void Loop_Force(param Param, particle* Particles){
  long i,j;
  double dx,dist;

  for (i=0; i<Param.N; i++){
    for (j=0; j<i; j++){
      dx = Particles[i].x - Particles[j].x;
      if (dx> Param.Lx_half) dx -= Param.Lx;
      if (dx<-Param.Lx_half) dx += Param.Lx;

      dist = dx;
      if (dx<0) dist=-dx;

      if (dx>1.0) continue;
      if (dx<-1.0) continue;

      Pair_force(dx, dist, i, j, Particles, Param);
    }
  }
}
#endif


#endif