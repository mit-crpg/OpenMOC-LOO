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
	double* _current_tot;
	double* _flux;
	double _normal;
	double _d_hat;
	double _d_tilde;
	int _mesh_cell;
	int _surface_num;
	int _crossings;
	double* _weight;
	int _id;

public:
	MeshSurface();
	virtual ~MeshSurface();
	void makeCurrents(int numAzim);
	meshSurfaceType getType();
	void setType(meshSurfaceType type);
	void setCurrent(double current, int group, int azim);
	double getCurrent(int group, int azim);
	void incrementCurrent(double current, int group, int azim);
	void setCurrentTot(double current, int group);
	double getCurrentTot(int group);
	void setNormal(double normal);
	double getNormal();
	void setDHat(double dHat);
	double getDHat();
	void setDTilde(double dTilde);
	double getDTilde();
	void setMeshCell(int meshCell);
	int getMeshCell();
	void setSurfaceNum(int surfaceNum);
	int getSurfaceNum();
	int getCrossings();
	void incrementCrossings(int crossings);
	void setCrossings(int crossings);
	double getWeight(int azim);
	void incrementWeight(double weight, int azim);
	void setWeight(double weight, int azim);
	int getId();
	void setId(int id);
};




#endif /* MESHSURFACE_H_ */
