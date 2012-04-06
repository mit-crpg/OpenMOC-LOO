/*
 * Plotter.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: samuelshaner
 */

#include "Plotter.h"


/**
 * Plotting constructor
 * @param geom a pointer to a geometry object
 */
Plotter::Plotter(Geometry* geom, const int bitDim, std::string extension, bool specs, bool fluxes) {

	/* extension for plotting files */
	_extension = extension;
	_geom = geom;
	_specs = specs;
	_fluxes = fluxes;

	_width = _geom->getWidth();
	_height = _geom->getHeight();
	double ratio = _width/_height;

	/* set pixel dimensions of plots */
	_bit_length_x = int (bitDim*ratio) + 1;
	_bit_length_y = bitDim + 1;

	/* make multiplier for translating geometry coordinates
	 * to bitmap coordinates.
	 */
	_x_pixel = double(_bit_length_x)/_width;
	_y_pixel = double(_bit_length_y)/_height;
}

/**
 * Plotter destructor frees memory for all tracks
 */
Plotter::~Plotter() {
}

/**
 * Return the x bit length of plots
 * @return the x bit length of plots
 */
int Plotter::getBitLengthX(){
	return _bit_length_x;
}

/**
 * Return the y bit length of plots
 * @return the y bit length of plots
 */
int Plotter::getBitLengthY(){
	return _bit_length_y;
}

/**
 * Return the x dimension pixel multiplier for plots
 * @return the x dimension pixel multiplier for plots
 */
double Plotter::getXPixel(){
	return _x_pixel;
}

/**
 * Return the y dimension pixel multiplier for plots
 * @return the y dimension pixel multiplier for plots
 */
double Plotter::getYPixel(){
	return _y_pixel;
}

/**
 * Return the plotting extension
 * @return the plotting extension
 */
std::string Plotter::getExtension(){
	return _extension;
}

/**
 * Return boolean to decide whether to plot specs
 * @return boolean to decide whether to plot specs
 */
bool Plotter::plotSpecs(){
	return _specs;
}

/**
 * Return boolean to decide whether to plot fluxes
 * @return boolean to decide whether to plot fluxes
 */
bool Plotter::plotFlux(){
	return _fluxes;
}

/* PLOTTING HELPER FUNCTIONS */
/* These functions manipulate data and call the generic plotting functions
 */

/**
 * Convert an x value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToBitmapX(double x){
	return int(x*_x_pixel + _bit_length_x/2);
}

/**
 * Convert an y value our from geometry coordinates to Bitmap coordinates.
 */
int Plotter::convertToBitmapY(double y){
	return int(-y*_y_pixel + _bit_length_y/2);
}

/**
 * Convert an x value our from geometry coordinates to Bitmap coordinates.
 */
double Plotter::convertToGeometryX(int x){
	return double((x - _bit_length_x/2.0) / _x_pixel);
}

/**
 * Convert an y value our from geometry coordinates to Bitmap coordinates.
 */
double Plotter::convertToGeometryY(int y){
	return double(-(y - _bit_length_y/2.0) / _y_pixel);
}

/**
 * Plots track segments in a pixMap array
 */
void Plotter::plotSegments(Track* track, double sin_phi,
		double cos_phi, int* pixMap){

	/* Initialize variables */
	double x0, y0, x1, y1;
	int num_segments;

	/* Set first segment start point and get the number of tracks*/
	x0 = track->getStart()->getX();
	y0 = track->getStart()->getY();
	num_segments = track->getNumSegments();

	log_printf(DEBUG, "looping over segments in plotter endy: %f...", track->getEnd()->getY());
	/* loop over segments and write to pixMap array */
	for (int k=0; k < num_segments; k++){
		x1 = x0 + cos_phi*track->getSegment(k)->_length;
		y1 = y0 + sin_phi*track->getSegment(k)->_length;

		/* "draw" segment on pixMap array */
		log_printf(DEBUG, "calling linefct x0: %f, y0: %f, x1: %f, y1: %f,length: %f, sin_phi: %f ...", x0, y0, x1, y1, track->getSegment(k)->_region_id, track->getSegment(k)->_length, sin_phi);
		LineFct(x0, y0, x1, y1, pixMap, track->getSegment(k)->_region_id);
		log_printf(DEBUG, "done calling linefct...");

		x0 = x1;
		y0 = y1;
	}
}

/**
 * plot given track and numReflect reflected tracks
 */
void Plotter::plotTracksReflective(Track* track, int numReflect){
	log_printf(NORMAL, "Writing tracks reflect bitmap...");

	/* initialize variables */
	double sin_phi, cos_phi, phi;
	Track *track2;
	bool get_out = TRUE;

	/* create BitMap for plotting */
	BitMap<int, double>* bitMap = new BitMap<int, double>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMap->color_type = RANDOM;
	bitMap->pixel_type = INTEGER;
	bitMap->pixels = new int[bitMap->pixel_x * bitMap->pixel_y];
	initialize(bitMap);

	/* loop through tracks and write to pixMap array */
	for (int i = 0; i < (numReflect + 1); i++){
		log_printf(DEBUG, "plotting reflective track: %d", i);
		log_printf(DEBUG, "x_start, y_start, x_end, y_end: %f, %f, %f, %f",
				track->getStart()->getX(),track->getStart()->getY(),
				track->getEnd()->getX(),track->getEnd()->getY());
		phi = track->getPhi();
		sin_phi = sin(phi);
		cos_phi = cos(phi);
		plotSegments(track, sin_phi, cos_phi, bitMap->pixels);

		/* Get next track */
		track2 = track;
		if (get_out){
			track = track->getTrackOut();
		}
		else {
			track = track->getTrackIn();
		}

		/*determine whether we want TrackIn or TrackOut of next track */
		if (track->getTrackOut() == track2){
			get_out = FALSE;
		}
		else{
			get_out = TRUE;
		}
	}

	/* plot pixMap array */
	plot(bitMap, "reflect", _extension);

	/* release memory */
	delete [] bitMap->pixels;
	delete bitMap;
}

/**
 * Bresenham's line drawing algorithm. Takes in the start and end coordinates
 * of line (in geometry coordinates), pointer to pixMap array, and line color.
 * "Draws" the line on pixMap array.
 * Taken from "Simplificaiton" code at link below
 * http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
*/
void Plotter::LineFct(double xIn, double yIn, double xOut, double yOut, int* pixMap, int color){

	/* initialize variables */
	int x0, y0, x1,y1;

	/* convert geometry coordinates to bitmap coordinates */
	x0 = convertToBitmapX(xIn);
	y0 = convertToBitmapY(yIn);
	x1 = convertToBitmapX(xOut);
	y1 = convertToBitmapY(yOut);

	log_printf(DEBUG, "linefct bitmap x0: %i, y0: %i, x1: %i, y1: %i", x0, y0, x1, y1);

	/* "draw" line on pixMap array */
	int dx = abs(x1-x0);
	int dy = abs(y1-y0);
	int sx, sy;
	if (x0 < x1){
		sx = 1;
	}
	else{
		sx = -1;
	}
	if (y0 < y1){
		sy = 1;
	}
	else{
		sy = -1;
	}
	int error = dx - dy;
	pixMap[y0 * _bit_length_x + x0] = color;
	pixMap[y1 * _bit_length_x + x1] = color;
	while (x0 != x1 && y0 != y1){
		pixMap[y0 * _bit_length_x + x0] = color;
		int e2 = 2 * error;
		if (e2 > -dy){
			error = error - dy;
			x0 = x0 + sx;
		}
		if (e2 < dx){
			error = error + dx;
			y0 = y0 + sy;
		}
	}
}

/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion) in either a
 * scaled Magick Plot or a silo plot.
 */
void Plotter::makeRegionMap(int* pixMapFSR, int* pixMap, int* regionMap){

	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = (int)regionMap[pixMapFSR[y * _bit_length_x + x]];
		}
	}
}


/**
 * Takes in a FSR pixMap array, a map (regionMap) that translates
 * a FSR id to a region id (cell, material, etc), and a description of the region
 * (regionName). Plots the resulting region pixMap (pixMapRegion) in either a
 * scaled Magick Plot or a silo plot.
 */
void Plotter::makeRegionMap(int* pixMapFSR, float* pixMap, double* regionMap){

	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			pixMap[y * _bit_length_x + x] = (float)regionMap[pixMapFSR[y * _bit_length_x + x]];
		}
	}
}


/**
 * Loops over pixels in pixMap array and finds the corresponding
 * FSR in geometry.
 */
void Plotter::makeFSRMap(int* pixMap){
	log_printf(NORMAL, "Generating FSR maps...");

	/* initialize variables */
	double x_global;
	double y_global;

	/* loop over pixels */
	for (int y=0;y< _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = convertToGeometryX(x);
			y_global = convertToGeometryY(y);

			log_printf(DEBUG, "finding cell for bit x: %i, bit y: %i, "
					"global x: %f, global %f", x, y, x_global, y_global);

			/* create point located in universe 0 */
			LocalCoords point(x_global,y_global);
			point.setUniverse(0);

			/* find which cell the point is in */
			_geom->findCell(&point);

			/* Store FSR id in pixMap */
			pixMap[y * _bit_length_x + x] = _geom->findFSRId(&point);

			/* Remove all allocated localcoords */
			point.prune();
		}
	}
}
