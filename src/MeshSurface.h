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


class MeshSurface {
private:
	double _current[NUM_ENERGY_GROUPS];
	double _flux[NUM_ENERGY_GROUPS];
	double _d_hat[NUM_ENERGY_GROUPS];
	double _d_tilde[NUM_ENERGY_GROUPS];
	double _d_dif[NUM_ENERGY_GROUPS];
	int _surface_num;
	int _id;
	int _cell_id;

public:
	MeshSurface();
	virtual ~MeshSurface();
	void makeCurrents();
	void setCurrent(double current, int group);
	double getCurrent(int group);
	void incrementCurrent(double current, int group);
	void setDHat(double dHat, int e);
	double* getDHat();
	void setDTilde(double dTilde, int e);
	double* getDTilde();
	void setDDif(double dTilde, int e);
	double* getDDif();
	void setSurfaceNum(int surfaceNum);
	int getSurfaceNum();
	int getId();
	void setId(int id);
	int getCellId();
	void setCellId(int id);

};


#endif /* MESHSURFACE_H_ */
