

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

