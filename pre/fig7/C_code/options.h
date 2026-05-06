
/* uncomment for non-interacting version */
// #define NI

/* random number generator */
#define _MT


/* threshold for converting to/from integer */
#define EPS 1e-10


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
#include "util/mt19937-64.c" // random number generator
#include "util/util.c" // utilities
#include "util/inhom_fluct.c" // inhomogeneous fluctuation landscape

#ifndef NI 
#include "util/spatial_hashing.c" // spatial hashing
#endif
