/*
 * MeshCell.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESHCELL_H_
#define MESHCELL_H_

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
#include "LocalCoords.h"
#include "MeshSurface.h"

class MeshCell {
private:
	double _width;
	double _height;
	std::vector<int> _FSRs;
	double _abs_rate;
	MeshSurface* _mesh_surfaces;
	double* _bounds;
	int _fsr_start;
	int _fsr_end;

public:
	MeshCell();
	virtual ~MeshCell();
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	std::vector<int>* getFSRs();
	void addFSR(int fsr);
	void setAbsRate(double absRate);
	double getAbsRate();
	MeshSurface* findSurface(LocalCoords* coord);
	void setBounds(double x, double y);
	double* getBounds();
	int getFSRStart();
	int getFSREnd();
	void setFSRStart(int fsr);
	void setFSREnd(int fsr);
	void makeSurfaces();
	MeshSurface* getMeshSurfaces(int surface);
};


#endif /* MESHCELL_H_ */
