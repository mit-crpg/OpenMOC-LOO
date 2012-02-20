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



class Plotter{
private:
	double _width;
	double _height;
	bool _plot_materials;
	bool _plot_cells;
	std::string _extension;
	Geometry* _geom;
	int _bit_length_x;
	int _bit_length_y;
	double _x_pixel;
	double _y_pixel;
	int* _pix_map_tracks;
	int* _pix_map_segments;
	int* _pix_map_FSR;
	int* _pix_map_materials;
	int* _pix_map_cells;
	std::map<int, std::string> _color_map;
public:
	Plotter(Geometry* geom, const int bit_dim, std::string extension,
			bool plotMaterials, bool plotCells);
	virtual ~Plotter();
	void plotMagick(int* pixMap, std::string type);
	void plotMagick(float* pixMap, std::string type);
	void plotTracksPng();
	void plotSegmentsBitMap(Track* track, double sin_phi, double cos_phi, int* map_array);
	void plotSegmentsPng();
	void LineFct(int x0, int y0, int x1, int y1, int* pixMap, int color = 1);
	void printPlottingTimers();
	void plotSilo(int* pixMap, std::string type);
	void plotSilo(float* pixMap, std::string type);
	void plotTracksReflective(Track* track, int numReflect);
	int* getPixMap(std::string type);
	void generateFsrMap();
	std::string getExtension();
	int getBitLengthX();
	int getBitLengthY();
	double getXPixel();
	double getYPixel();
	void plot(int* pixMap, std::string type);
	void plot(float* pixMap, std::string type);
};





#endif /* PLOTTER_H_ */
