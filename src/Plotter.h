/*
 * Plotter.h
 *
 *  Created on: Feb 15, 2012
 *      Author: samuelshaner
 */

#ifndef PLOTTER_H_
#define PLOTTER_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "Point.h"
#include "Track.h"
#include "Geometry.h"
#include "Magick++.h"
#include "silo.h"
#include "LocalCoords.h"
#include "Cell.h"
#include "quickplot.h"
#include "Mesh.h"
#include "MeshCell.h"

class Plotter{
private:
	double _width;
	double _height;
	std::string _extension;
	Geometry* _geom;
	int _bit_length_x;
	int _bit_length_y;
	double _x_pixel;
	double _y_pixel;
	bool _specs;
	bool _fluxes;
	bool _net_current;
	bool _plot_diffusion;
	bool _plot_keff;
public:
	Plotter(Geometry* geom, const int bit_dim, std::string extension, bool specs, bool fluxes,
			bool netCurrent, bool plotDiffusion, bool plotKeff);
	virtual ~Plotter();
	void plotTracksReflective(Track* track, int numReflect);
	void makeFSRMap(int* pixMap);
	int getBitLengthX();
	int getBitLengthY();
	double getXPixel();
	double getYPixel();
	bool plotSpecs();
	bool plotFlux();
	std::string getExtension();
	void makeRegionMap(int* pixMapFSR, int* pixMap, int* regionMap);
	void makeRegionMap(int* pixMapFSR, float* pixMap, double* regionMap);
	double convertToGeometryX(int x);
	double convertToGeometryY(int y);
	void plotCMFDMesh(Mesh* mesh);
	int convertToPixelX(double x);
	int convertToPixelY(double y);
	void plotNetCurrents(Mesh* mesh);
//	void plotSurfaceFlux(Mesh* mesh);
	void plotDHats(Mesh* mesh, int iter_num);
	void plotXS(Mesh* mesh, int iter_num);
	bool plotCurrent();
	bool plotKeff();
	bool plotDiffusion();
	void plotCMFDflux(Mesh* mesh, std::string string, int iter_num);
	void plotCMFDKeff(Mesh* mesh, int num_iter);
};





#endif /* PLOTTER_H_ */
