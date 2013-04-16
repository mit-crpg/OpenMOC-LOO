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
#include "Surface.h"



class MeshSurface {
private:
	double _current[NUM_ENERGY_GROUPS];
	double _flux[NUM_ENERGY_GROUPS][2];
	double _d_hat[NUM_ENERGY_GROUPS];
	double _d_tilde[NUM_ENERGY_GROUPS];
	double _d_dif[NUM_ENERGY_GROUPS];

	/* Surface ID could be 0,1,2,3 and is later set in MeshCell.cpp */
	int _id; 
	int _cell_id;
	boundaryType _boundary_type;

public:
	MeshSurface();
	virtual ~MeshSurface();

	/* LOO Only */
	void setFlux(double flux, int group, int index);
	double getFlux(int group, int index);
	void incrementFlux(double flux, int group, int index);

	/* CMFD Only */
	void makeCurrents();
	void setCurrent(double current, int group);
	double getCurrent(int group);
	void incrementCurrent(double *current);
	void setDHat(double dHat, int e);
	double* getDHat();
	void setDTilde(double dTilde, int e);
	double* getDTilde();
	void setDDif(double dTilde, int e);
	double* getDDif();

	/* General Purpose */
	int getId();
	void setId(int id);
	int getCellId();
	void setCellId(int id);
	void setBoundary(boundaryType boundary);
	boundaryType getBoundary();
};


#endif /* MESHSURFACE_H_ */
