
/* spatial dimension */
#define TWOD


/* whether to use temperature landscape */
#define TOFX_CUB
// #define TOFX_LIN
// #define TOFX_TRIG


/* whether to use force landscape */
// #define FOFX_CUB
// #define FOFX_LIN

/* what, if any, force to use */
//#define FORCE_LINEAR
#define FORCE_HARMONIC 
//#define FORCE_IHARMONIC
// #define FORCE_CUBIC
//#define FORCE_QUARTIC
//#define FORCE_WCA 


/* what data storage routines to do */
#define STORE_DISP
#define STORE_PROF
// #define STORE_JPROF
#define STORE_SIGMAIKPROF
#define STORE_TPROF
// #define STORE_FPROF
//#define TRACK
#define STORE_EPR


/*
  Boundary condition (periodic or closed)

  If closed boundary conditions are used, then a confining potential
  is inserted INSIDE the box, to keep the particles far away from the
  periodic boundary condition of the cell. If particles get close
  enough to the wall, the program protests.
*/
//#define CBC_X
//#define CBC_Y


/* initial condition: random, given, or slab */
#define RANDOM_IC
//#define GIVEN_IC
//#define SLAB_IC


/* random number generator */
#define _MT
//#define _PCG


/* threshold for converting to/from integer */
#define EPS 1e-10


////////////////////////////////////////////////////////////////////////////////

/* 
  Options that are consequences of other options 
*/

#if defined(FOFX_LIN) || defined(FOFX_CUB) || defined(FOFX_TRIG)
#define FOFX
#endif

#ifdef FOFX_CUB
#define FLUCT_CUB
#elif defined(FOFX_LIN)
#define FLUCT_LIN
#endif

#if defined(TOFX_LIN) || defined(TOFX_CUB) || defined(TOFX_TRIG)
#define TOFX
#endif

#ifdef TOFX_CUB
#define FLUCT_CUB
#elif defined(TOFX_LIN)
#define FLUCT_LIN
#elif defined(TOFX_TRIG)
#define FLUCT_TRIG
#endif

#if !defined(FORCE_LINEAR) & !defined(FORCE_HARMONIC) & !defined(FORCE_IHARMONIC) & !defined(FORCE_CUBIC) & !defined(FORCE_QUARTIC) & !defined(FORCE_WCA)
#define NI
#endif

#if defined(STORE_PROF) || defined(STORE_JPROF) || defined(STORE_SIGMAIKPROF) || defined(STORE_FPROF) || defined(STORE_TPROF)
#define STORE_ANYPROF
#endif


////////////////////////////////////////////////////////////////////////////////

/*
  Include all dependencies
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // Library that measures dates and times
#include <assert.h>
#include "mt19937-64.c" // random number generator
#include "util.c" // utilities

#ifndef NI 
#include "spatial_hashing.c" // spatial hashing
#endif

#if defined(TOFX) || defined(FOFX)
#include "inhom_fluct.c" // inhomogeneous fluctuation landscape
#endif
