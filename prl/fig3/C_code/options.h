
/* random number generator */
#define _MT

/* threshold for converting to/from integer */
#define EPS 1e-10

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
#include "spatial_hashing.c" // spatial hashing
#include "inhom_fluct.c" // inhomogeneous fluctuation landscape