/*
 * configurations.h
 *
 *  Created on: Jan 16, 2012
 *      Author: will
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_

#include <math.h>

/* Simulation specific */
#define NUM_POLAR_ANGLES 3
#define NUM_ENERGY_GROUPS 7
#define GRP_TIMES_ANG NUM_POLAR_ANGLES*NUM_ENERGY_GROUPS

/* Constants */
#define FOUR_PI 12.5663706143
#define ONE_OVER_FOUR_PI 0.0795774715

/* Error thresholds */
#define ON_SURFACE_POS 1E-3
#define ON_SURFACE_NEG -1E-3

#endif /* CONFIGURATIONS_H_ */
