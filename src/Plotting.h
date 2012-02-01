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
	int _x_pixel;
	int _y_pixel;
	Geometry* _geom;
	std::list<Magick::Drawable> _draw_lists[5];
public:
	Plotting(Geometry* geom);
	virtual ~Plotting();
	void plotTracksTiff(TrackGenerator* track_generator);
	// void addTrackSegments(Track* track);
	void plotSegments(TrackGenerator* track_generator);
};


#endif /* PLOTTING_H_ */
