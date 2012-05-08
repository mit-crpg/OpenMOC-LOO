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
#include "log.h"


class Mesh {
private:
	MeshCell* _cells;
	int _cell_width;
	int _cell_height;
	double _width;
	double _height;

public:
	Mesh();
	virtual ~Mesh();
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	int getCellWidth();
	int getCellHeight();
	void setCellWidth(int cellWidth);
	void setCellHeight(int cellHeight);
	void makeMeshCells();
	MeshCell* getCells();
	int findMeshCell(double x, double y);

};

#endif /* MESH_H_ */
