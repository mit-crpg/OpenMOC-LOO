/*
 * Track.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Track.h"


/**
 * Track constructor
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 */
Track::Track(double start_x, double start_y, double end_x,
		double end_y, double phi) {
	_start.setCoords(start_x, start_y);
	_end.setCoords(end_x, end_y);
	_phi = phi;
}


Point Track::getEnd() const {
    return _end;
}

double Track::getPhi() const
{
    return _phi;
}

Point Track::getStart() const
{
    return _start;
}

double Track::getWeight() const
{
    return _weight;
}

/**
 * Track destructor
 */
Track::~Track() { }




