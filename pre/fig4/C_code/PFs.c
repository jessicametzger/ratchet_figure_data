

/*
  HARMONIC: We use harmonic spheres F = -\grad V where
  
  V(r)= (amp sigma / 2) (1 - r/sigma)^2

  and force magnitude

  F(r_i-r_j) = amp (1 - |r_i-r_j|/sigma)
  
  We use a cut-off at r=sigma 
*/
double get_force_harmonic_params(double k, double a, double mu, double dt, double* c1, double* c2){
  c1[0] = k * dt * mu;
  c2[0] = -k * dt * mu / a;
}
double Force_harmonic(double dist, double c1, double c2){
	return c1 + c2*dist;
}

// given (corrected) distances dx,dy and indices j,k, update forces (stored in "Particles").
// ASSUME THEY ARE CLOSE ENOUGH TO INTERACT!
void Pair_force(double dx, double dy, double dist2, long j, long k, particle* Particles, param Param){
  double Force,Forcex,Forcey;

  double dist=pow(dist2,0.5);

  Force = Force_harmonic(dist,Param.c1,Param.c2);
  Forcex = Force*dx/dist;
  Forcey = Force*dy/dist;
  Particles[j].ifxs[k] = Forcex;
  Particles[j].ifys[k] = Forcey;
  Particles[k].ifxs[j] = -Forcex;
  Particles[k].ifys[j] = -Forcey;
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

    x = Particles[j].x - Particles[j].dx;
    y = Particles[j].y - Particles[j].dy;
    while(x>Param.Lx) x-=Param.Lx;
    while(x<0) x+=Param.Lx;
    while(y>Param.Ly) y-=Param.Ly;
    while(y<0) y+=Param.Ly;
    bi = (int) (floor(x/Param.rmax) + EPS);
    bj = (int) (floor(y/Param.rmax) + EPS);
    
    box_left = Param.rmax * floor(x/Param.rmax);
    box_bottom = Param.rmax * floor(y/Param.rmax);
    box_right = box_left + Param.rmax;
    box_top = box_bottom + Param.rmax;

    for (k=0; k<Param.N; k++){
      if (k==j) continue;
      xn = Particles[k].x - Particles[k].dx;
      yn = Particles[k].y - Particles[k].dy;
      while(xn>Param.Lx) xn-=Param.Lx;
      while(xn<0) xn+=Param.Lx;
      while(yn>Param.Ly) yn-=Param.Ly;
      while(yn<0) yn+=Param.Ly;
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
