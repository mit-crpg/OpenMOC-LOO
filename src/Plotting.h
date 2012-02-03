/*
 * Plotting.h
 *
 *  Created on: Feb 1, 2012
 *      Author: samuelshaner
 */

#ifndef PLOTTING_H_
#define PLOTTING_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include "Point.h"
#include "Track.h"
#include "Geometry.h"
#include "TrackGenerator.h"
#include "Magick++.h"

class TrackGenerator;

class Plotting {
private:
	int _bit_length_x;
	int _bit_length_y;
	double _x_pixel;
	double _y_pixel;
	int* _pix_map_tracks;
	int* _pix_map_segments;
	Geometry* _geom;
	std::list<Magick::Drawable> _draw_lists[7];
public:
	Plotting(Geometry* geom, int bitDim);
	virtual ~Plotting();
	void plotTracksTiff(TrackGenerator* track_generator);
	// void addTrackSegments(Track* track);
	void plotSegments(TrackGenerator* track_generator);
	int sgn (long a);
	void LineFct(int a, int b, int c, int d, int col);
	void SegFct(int a, int b, int c, int d, int col);

};


#endif /* PLOTTING_H_ */
