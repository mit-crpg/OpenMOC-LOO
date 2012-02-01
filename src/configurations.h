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
#define NUM_ENERGY_GROUPS 7
#define GRP_TIMES_ANG NUM_POLAR_ANGLES*NUM_ENERGY_GROUPS



/******************************************************************************
 *********************** PHYSICAL CONSTANTS ***********************************
 *****************************************************************************/

#define FOUR_PI 12.5663706143
#define ONE_OVER_FOUR_PI 0.0795774715



/******************************************************************************
 *************************** ERROR THRESHOLDS *********************************
 *****************************************************************************/

#define ON_SURFACE_THRESH 1E-8
#define ON_LATTICE_CELL_THRESH 1E-6
#define TINY_MOVE 1E-6



#endif /* CONFIGURATIONS_H_ */
