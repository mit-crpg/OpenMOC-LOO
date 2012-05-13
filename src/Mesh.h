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
#include "LocalCoords.h"
#include "MeshSurface.h"


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
	void setCellBounds();
	void setFSRBounds();
	int findMeshCell(double x, double y);
	MeshSurface* findMeshSurface(int fsr_id, LocalCoords* coord);
	void printBounds();

};

#endif /* MESH_H_ */
