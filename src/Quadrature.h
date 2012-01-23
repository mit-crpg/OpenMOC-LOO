/*
 * Quadrature.h
 *
 *  Created on: Jan 20, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <stdlib.h>
#include <iostream>
#include "configurations.h"
#include "log.h"

enum quadratureType {
	LEONARD,
	TABUCHI
};

class Quadrature {
private:
	quadratureType _type;
	double _sinthetas[NUM_POLAR_ANGLES];
	double _weights[NUM_POLAR_ANGLES];
	double _multiples[NUM_POLAR_ANGLES];
public:
	Quadrature(quadratureType type);
	virtual ~Quadrature();
	quadratureType getType();
	double getSinTheta(int n);
	double getWeight(int n);
	double getMultiple(int n);
	double* getSinThetas();
    double* getWeights();
    double* getMultiples();
};

#endif /* QUADRATURE_H_ */
