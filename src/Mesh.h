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
#include "configurations.h"


class Mesh {
private:
	MeshCell* _cells;
	int _cell_width;
	int _cell_height;
	double _width;
	double _height;
	bool _multigroup;
	bool _print_matrices;
	boundaryType _boundary[4];
	double _keff_cmfd[int(MAX_ITERATIONS)];
	double _keff_moc[int(MAX_ITERATIONS)];
	double _old_time;
	int* _fsr_indices;
	double* _cell_bounds;
	MeshSurface **_surfaces;

public:
	Mesh();
	virtual ~Mesh();
	void makeMeshCells();
	double getWidth();
	double getHeight();
	void setWidth(double width);
	void setHeight(double height);
	int getCellWidth();
	int getCellHeight();
	void setCellWidth(int cellWidth);
	void setCellHeight(int cellHeight);
	MeshCell* getCells(int cell);
	void setCellBounds();
	void setFSRBounds(boundaryType left, boundaryType right, boundaryType bottom, boundaryType top);
	int findMeshCell(double x, double y);
	int findMeshSurface(int fsr_id, LocalCoords* coord);
	void printBounds();
	void printCurrents();
	void computeTotCurrents();
	void splitCorners();
	void setMultigroup(bool multigroup);
	bool getMultigroup();
	void setPrintMatrices(bool printMatrices);
	bool getPrintMatrices();
	void setBoundary(boundaryType boundary, int s);
	boundaryType getBoundary(int s);
	void setKeffCMFD(double keff, int iter);
	double getKeffCMFD(int iter);
	void setKeffMOC(double keff, int iter);
	double getKeffMOC(int iter);
	MeshSurface **getSurfaces();
	double getOldTime();
	void setOldTime(double time);

};

#endif /* MESH_H_ */
