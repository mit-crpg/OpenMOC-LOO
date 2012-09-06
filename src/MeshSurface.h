/*
 * MeshSurface.h
 *
 *  Created on: May 12, 2012
 *      Author: samuelshaner
 */

#ifndef MESHSURFACE_H_
#define MESHSURFACE_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <list>
#include "log.h"
#include <vector>
#include "configurations.h"

enum meshSurfaceType{
	SIDEX,
	SIDEY,
	CORNER
};

class MeshSurface {
private:
	meshSurfaceType _type;
	double* _current;
	double* _flux;
	double _normal;
	double _d_hat;
	double _d_tilde;

public:
	MeshSurface();
	virtual ~MeshSurface();
	meshSurfaceType getType();
	void setType(meshSurfaceType type);
	void setCurrent(double current, int group);
	double getCurrent(int group);
	void incrementCurrent(double current, double phi, int group);
	void setNormal(double normal);
	double getNormal();
	void setFlux(double flux, int group);
	double getFlux(int group);
	void incrementFlux(double flux, int group);
	void setDHat(double dHat);
	double getDHat();
	void setDTilde(double dTilde);
	double getDTilde();
};



#endif /* MESHSURFACE_H_ */
