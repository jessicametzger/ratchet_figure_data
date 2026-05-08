#define _MT
//#define _PCG

#define EPS 1e-10


//#define TRACK


/*
  Choose the initial condition (IC) between random IC, if RANDOMIC is
  defined, or read from an input if GIVENIC is defined
*/
#define RANDOM_IC


/*
  Choose the boundary conditions between periodic, if PBC is defined,
  and closed along x, if CLOSEDBC is defined
  
  If closed boundary conditions are used, then a confining potential
  is inserted INSIDE the box, to keep the particles far away from the
  periodic boundary condition of the cell. If particles get close
  
  enough to the wall, the program protests.
*/
#define PBC


