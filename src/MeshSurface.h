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
#include <map>

enum meshSurfaceType{
	SIDE,
	CORNER
};

class MeshSurface {
private:
	meshSurfaceType _type;
	double _current;

public:
	MeshSurface();
	virtual ~MeshSurface();
	meshSurfaceType getType();
	void setType(meshSurfaceType type);
};


#endif /* MESHSURFACE_H_ */
