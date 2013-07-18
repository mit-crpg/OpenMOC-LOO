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
#include <sstream>
#include <string>
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
    quadratureType getType() const;
    double getSinTheta(const int n) const;
    double getWeight(const int n) const;
    double getMultiple(const int n) const;
    double* getSinThetas();
    double* getWeights();
    double* getMultiples();
    std::string toString();
};

#endif /* QUADRATURE_H_ */
