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

public:
	MeshSurface();
	virtual ~MeshSurface();
	meshSurfaceType getType();
	void setType(meshSurfaceType type);
	void setCurrent(double current, int group);
	double getCurrent(int group);
	void incrementCurrent(double current, int group);
};


#endif /* MESHSURFACE_H_ */
