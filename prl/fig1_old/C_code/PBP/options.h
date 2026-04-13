#define _MT
//#define _PCG

#define EPS 1e-10

//#define TRACK

// #define EPR

#define SIGMA


#define NI


//#define FORCE_LINEAR
/*
  "tent" force
  
  V(r) = amp sigma (1-|r|/sigma)
  
  F(r) = amp 
  
  We use a cut-off at r=sigma
*/

//#define FORCE_QUARTIC
/*
  Quartic force F = -grad V where

  V(r) = amp*(1 + x^4/sigma^4 - 2 x^2/sigma^2)
*/

//#define HARMONIC 
/*
  We use V(r) = (amp sigma / 2) (1 - (r/sigma)^2)
  
  and force magnitude
  
  F(r_i-r_j) = amp |r_i-r_j| / sigma
  
  We use a cut-off at r=sigma
  
*/

// #define FORCE_INV_HARMONIC
/*
  We use harmonic spheres F = -\grad V where
  
  V(r)= (amp sigma / 2) (1 - (r/sigma)^2)

  and force magnitude

  F(r_i-r_j) = amp (1 - |r_i-r_j|/sigma)
  
  We use a cut-off at r=sigma 
*/

//#define WCA 
/*
  We use WCA potential for F = -\grad V where
  
  V=4 epsilon (sigma^12/r^12-sigma^6/r^6) + epsilon
  
  so the force magnitude is
  
  F = 4 epsilon (12 sigma^12 / r^13 - 6 sigma^6 / r^7)
  
  We use the LJ potential with a cut-off at r=sigma*2^1/6
  
  There is for now no adaptative time-stepping because of the order dt and \sqrt{dt}.
*/



/*
  Choose the initial condition (IC) between random IC, if RANDOMIC is
  defined, or read from an input if GIVENIC is defined
  TO ADD: SLAB_IC
*/
#define RANDOM_IC
//#define GIVEN_IC


/*
  Choose the boundary conditions between periodic, if PBC is defined,
  and closed along x, if CLOSEDBC is defined
  
  If closed boundary conditions are used, then a confining potential
  is inserted INSIDE the box, to keep the particles far away from the
  periodic boundary condition of the cell. If particles get close
  enough to the wall, the program protests.
  
  TO ADD: CLOSEDBC
*/
#define PBC

/* 
  whether or not we want x-dependent temperature, mobility, 
  and external potential fields.
*/
#define TOFX
//#define MUOFX
//#define UOFX

// use non-Ito discretizatiom
// #define ALPHA_NEQ_0

/*
  what type of landscape for T(x), mu(x), U(x).
  Linear is linear interpolation between two platforms.
  Cubic is cubic interpolation between two platforms (differentiable everywhere).
  Trig is sum of two trig functions in x,y variables.
*/
//#define LANDSCAPE_CUBIC
// #define LANDSCAPE_LINEAR
#define LANDSCAPE_TRIG


#define OVERDAMP
//#define UNDERDAMP


/*
  Which integrator we will use
*/
#define EULER

