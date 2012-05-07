/*
 * Mesh.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESH_H_
#define MESH_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "MeshCell.h"


class Mesh {
private:
	MeshCell* _cells;
	int _width;
	int _height;

public:
	Mesh(int width, int height);
	virtual ~Mesh();
	int getWidth();
	int getHeight();
	void setWidth(int width);
	void setHeight(int height);
	MeshCell* getCells();

};

#endif /* MESH_H_ */
