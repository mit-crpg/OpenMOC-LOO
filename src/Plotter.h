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
#include "Point.h"
#include "Track.h"
#include "Geometry.h"
#include "Magick++.h"
#include "silo.h"
#include "LocalCoords.h"
#include "Cell.h"
#include <visit_writer.h>



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
	std::map<int, std::string> _color_map;
public:
	Plotter(Geometry* geom, const int bit_dim, std::string extension);
	virtual ~Plotter();
	void plotMagick(int* pixMap, std::string type);
	void plotMagick(float* pixMap, std::string type);
	void plotSegments(Track* track, double sin_phi, double cos_phi, int* pixMap);
	void LineFct(double x0, double y0, double x1, double y1, int* pixMap, int color = 1);
	void plotSilo(int* pixMap, std::string type);
	void plotSilo(float* pixMap, std::string type);
	void plotTracksReflective(Track* track, int numReflect);
	void plotFSRs(int* pixMap);
	int getBitLengthX();
	int getBitLengthY();
	double getXPixel();
	double getYPixel();
	void plot(int* pixMap, std::string type);
	void plot(float* pixMap, std::string type);
	void FlipBitmap(int* pixMap);
	void FlipBitmap(float* pixMap);
	int convertToBitmapX(double x);
	int convertToBitmapY(double y);
	void plotRegion(int* pixMap, int* regionMap, std::string regionName);
	void plotRegion(int* pixMap, double* regionMap, std::string regionName);
	void initializePixMap(int* pixMap);
	void initializePixMap(float* pixMap);
	double convertToGeometryX(int x);
	double convertToGeometryY(int y);
};





#endif /* PLOTTER_H_ */
