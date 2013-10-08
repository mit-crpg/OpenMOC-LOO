/*
 * configurations.h
 *
 *  Created on: Jan 16, 2012
 *      Author: will
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_

#include <math.h>

/****** For debugging, NEW means the new implementation with weights *****/
#define NEW 1



/******************************************************************************
 ****************************** USER DEFINED **********************************
 *****************************************************************************/
#define NUM_POLAR_ANGLES 1
#define P0 0.798184
#define NUM_ENERGY_GROUPS 7
#define GRP_TIMES_ANG NUM_POLAR_ANGLES*NUM_ENERGY_GROUPS

/* The number of old k_eff values which will be stored. Keeping track of more
 * than just the most recent k_eff is intended to prevent from prematurely
 * exiting fixed source iteration at a local minimum or maximum */
#define NUM_KEFFS_TRACKED 3

/** Maximum number of fixed source iterations allowed */
#define MAX_ITERATIONS 2e3

/* Precompute and store exponential pre-factors in transport equation */
#define STORE_PREFACTORS true

/* Number of significant digits for computing hashmap exponential prefactors */
#define FSR_HASHMAP_PRECISION 12

/* If this machine has OpenMP installed, define as true for parallel speedup */
#define USE_OPENMP false

/******************************************************************************
 *********************** PHYSICAL CONSTANTS ***********************************
 *****************************************************************************/
#define FOUR_PI 12.566370614359172
#define TWO_PI 6.283185307179586
#define ONE_OVER_FOUR_PI 0.07957747154594767
#define PI 3.141592653589793
#define SIN_THETA_45 0.70710678118654746

/******************************************************************************
 *************************** ERROR THRESHOLDS *********************************
 *****************************************************************************/
/* Error threshold for determining how close a point needs to be to a surface
 * to be considered on it */
#define ON_SURFACE_THRESH 1E-12

/* Error threshold for determining how close to the boundary of a lattice cell
 * a point needs to be to be considered on it */
#define ON_LATTICE_CELL_THRESH 1E-12

/* Distance a point is moved to cross over a surface into a new cell during
 * track segmentation */
#define TINY_MOVE 1E-10

#define TINY_ANGLE 1E-10

#endif /* CONFIGURATIONS_H_ */
