/*
 * configurations.h
 *
 *  Created on: Jan 16, 2012
 *      Author: will
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_

#include <math.h>

/******************************************************************************
 ****************************** USER DEFINED **********************************
 *****************************************************************************/

#define NUM_POLAR_ANGLES 3
#define NUM_ENERGY_GROUPS 2
#define GRP_TIMES_ANG NUM_POLAR_ANGLES*NUM_ENERGY_GROUPS

/* Convergence threshold for scalar flux in each region during fixed source
 * iteration */
#define FLUX_CONVERGENCE_THRESH 1E-6

/* The number of old k_eff values which will be stored. Keeping track of more
 * than just the most recent k_eff is intended to prevent from prematurely
 * exiting fixed source iteration at a local minimum or maximum */
#define NUM_KEFFS_TRACKED 3

/** Maximum number of fixed source iterations allowed */
#define MAX_ITERATIONS 3E+3

/* Precompute and store exponential pre-factors in transport equation */
#define STORE_PREFACTORS false

/* Number of significant digits for computing hashmap exponential prefactors */
#define FSR_HASHMAP_PRECISION 6

/* If this machine has OpenMP installed, define as true for parallel speedup */
#define USE_OPENMP false

/* Axial buckling term used for diffusion problems */
#define AXIAL_BUCKLING 1e-4

/******************************************************************************
 *********************** PHYSICAL CONSTANTS ***********************************
 *****************************************************************************/

#define FOUR_PI 12.5663706143
#define ONE_OVER_FOUR_PI 0.0795774715
#define PI 3.1415926536


/******************************************************************************
 *************************** ERROR THRESHOLDS *********************************
 *****************************************************************************/

/* Error threshold for determining how close the sum of sigma_a and sigma_s
 * must match that of sigma_t for each energy group */
#define SIGMA_T_THRESH 1E-3

/* Error threshold for determining how close a point needs to be to a surface
 * to be considered on it */
#define ON_SURFACE_THRESH 1E-12

/* Error threshold for determining how close to the boundary of a lattice cell
 * a point needs to be to be considered on it */
#define ON_LATTICE_CELL_THRESH 1E-12

/* Distance a point is moved to cross over a surface into a new cell during
 * track segmentation */
#define TINY_MOVE 1E-10


#endif /* CONFIGURATIONS_H_ */
