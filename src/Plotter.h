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
#include "plotterNew.h"



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
public:
	Plotter(Geometry* geom, const int bit_dim, std::string extension, bool specs, bool fluxes);
	virtual ~Plotter();
	void plotSegments(Track* track, double sin_phi, double cos_phi, int* pixMap);
	void LineFct(double x0, double y0, double x1, double y1, int* pixMap, int color = 1);
	void plotTracksReflective(Track* track, int numReflect);
	void makeFSRMap(int* pixMap);
	int getBitLengthX();
	int getBitLengthY();
	double getXPixel();
	double getYPixel();
	bool plotSpecs();
	bool plotFlux();
	std::string getExtension();
	int convertToBitmapX(double x);
	int convertToBitmapY(double y);
	void makeRegionMap(int* pixMapFSR, int* pixMap, int* regionMap);
	void makeRegionMap(int* pixMapFSR, float* pixMap, double* regionMap);
	double convertToGeometryX(int x);
	double convertToGeometryY(int y);
};





#endif /* PLOTTER_H_ */
