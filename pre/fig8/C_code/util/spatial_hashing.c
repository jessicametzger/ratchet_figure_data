#ifndef ONED

// Structure used to compute forces in the presence of periodic boundary conditions
typedef struct box{
  long i; // position of box along x axis
  long j; // position of box along y axis
  double epsilonx; //distance to add along x to take into account PBC
  double epsilony; //distance to add along y to take into account PBC
} box;

void AddinBox(long i,long bi, long bj,long** Boxes,long* Neighbours){
  long k;
  // Save the old state of Box, -1 or a particle number
  k=Boxes[bi][bj];
  
  //Put particle i at the top
  Boxes[bi][bj]=i;
  // There is no one above particle i
  Neighbours[2*i]=-1;
  // k is after i
  Neighbours[2*i+1]=k;
  // If k is particle, i is before k
  if(k!=-1)
    Neighbours[2*k]=i;
}

/* 
   To remove someone from the box, you have to remove it, but you do
   not need to update its neighbours. This will be done when it is
   added to a new box
*/
void RemovefromBox(long i,long bi,long bj,long** Boxes,long* Neighbours){
  long next;
  long prev;
  // Store the next particle
  next=Neighbours[2*i+1];
  // If the particle is the first in the box
  if(Boxes[bi][bj]==i){
    // next is the new top of the box
    Boxes[bi][bj]=next;
    // if next is a particle, there is no one before it
    if(next!=-1){
      Neighbours[2*next]=-1;
    }
  }
  // If there is someone before you
  else{
    // Store the previous particle
    prev=Neighbours[2*i];
    // The next of the previous is your next
    Neighbours[2*prev+1]=next;
    // If next is a particle, its previous is your previous
    if(next!=-1)
      Neighbours[2*next]=prev;
  }
}


void DefineNeighbouringBoxes(box *** NeighbouringBoxes, int NxBox, int NyBox, double Lx, double Ly){
  int i,j;
  for (i=0;i<NxBox;i++){
    NeighbouringBoxes[i] = (box**) malloc(NyBox*sizeof(box*));
    for (j=0;j<NyBox;j++){
      NeighbouringBoxes[i][j] = (box*) malloc(5*sizeof(box));
      
      // This is the same box
      NeighbouringBoxes[i][j][0].i=i;
      NeighbouringBoxes[i][j][0].j=j;
      NeighbouringBoxes[i][j][0].epsilonx=0;
      NeighbouringBoxes[i][j][0].epsilony=0;
      
      // The box above
      NeighbouringBoxes[i][j][1].i=i;
      NeighbouringBoxes[i][j][1].j=(j+1)%NyBox;
      NeighbouringBoxes[i][j][1].epsilonx=0;
      NeighbouringBoxes[i][j][1].epsilony=(j==NyBox-1)?(Ly):0;
      
      // The box above, to the right
      NeighbouringBoxes[i][j][2].i=(i+1)%NxBox;
      NeighbouringBoxes[i][j][2].j=(j+1)%NyBox;
      NeighbouringBoxes[i][j][2].epsilonx=(i==NxBox-1)?(Lx):0;
      NeighbouringBoxes[i][j][2].epsilony=(j==NyBox-1)?(Ly):0;
      
      // The box to the right
      NeighbouringBoxes[i][j][3].i=(i+1)%NxBox;
      NeighbouringBoxes[i][j][3].j=j;
      NeighbouringBoxes[i][j][3].epsilonx=(i==NxBox-1)?(Lx):0;
      NeighbouringBoxes[i][j][3].epsilony=0;
      
      // The box below, to the right
      NeighbouringBoxes[i][j][4].i=(i+1)%NxBox;
      NeighbouringBoxes[i][j][4].j=(j-1+NyBox)%NyBox;
      NeighbouringBoxes[i][j][4].epsilonx=(i==NxBox-1)?(Lx):0;
      NeighbouringBoxes[i][j][4].epsilony=(j==0)?(-Ly):0;
    }
  }
}




#else

// Structure used to compute forces in the presence of periodic boundary conditions
typedef struct box{
  long i; // position of box along x axis
  double epsilonx; //distance to add along x to take into account PBC
} box;

/*
  To add someone in a box, you need to put it there and to set its new
  neighbours. This functions add at the top.
*/
void AddinBox(long i,long bi, long* Box,long* Neighbours){
  long k;
  // Save the old state of Box, -1 or a particle number
  k=Box[bi];
  
  //Put particle i at the top
  Box[bi]=i;
  // There is no one above particle i
  Neighbours[2*i]=-1;
  // k is after i
  Neighbours[2*i+1]=k;
  // If k is particle, i is before k
  if(k!=-1)
    Neighbours[2*k]=i;
}


/* 
   To remove someone from the box, you have to remove it, but you do
   not need to update its neighbours. This will be done when it is
   added to a new box
*/
void RemovefromBox(long i,long bi,long* Box,long* Neighbours){
  long next;
  long prev;
  // Store the next particle
  next=Neighbours[2*i+1];
  // If the particle is the first in the box
  if(Box[bi]==i){
    // next is the new top of the box
    Box[bi]=next;
    // if next is a particle, there is no one before it
    if(next!=-1){
      Neighbours[2*next]=-1;
    }
  }
  // If there is someone before you
  else{
    // Store the previous particle
    prev=Neighbours[2*i];
    // The next of the previous is your next
    Neighbours[2*prev+1]=next;
    // If next is a particle, its previous is your previous
    if(next!=-1)
      Neighbours[2*next]=prev;
  }
}



void DefineNeighbouringBoxes(box ** NeighbouringBoxes,int NxBox, double Lx){
  int i;
  for (i=0;i<NxBox;i++){
    NeighbouringBoxes[i] = (box*) malloc(2*sizeof(box));
    
    NeighbouringBoxes[i][0].i = i;
    NeighbouringBoxes[i][0].epsilonx = 0;

    NeighbouringBoxes[i][1].i = (i+1)%NxBox;
    NeighbouringBoxes[i][1].epsilonx = (i==NxBox-1)?(Lx):0;
  }
}
#endif