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


/**
 * Track destructor
 */
Track::~Track() {
	clearSegments();
}


/**
 * Set the track azimuthal weight
 * @param weight the azimuthal weight
 */
void Track::setWeight(double weight) {
    _weight = weight;
}


/**
 * Adds a segment pointer to this Track's list of segments
 * IMPORTANT: assumes that segments are added in order of their starting
 * location from the track's start point
 * @param segment a pointer to the segment
 */
void Track::addSegment(segment* segment) {
	_segments.push_back(segment);
}


/**
 * Returns the track's end point
 * @return a pointer to the track's end point
 */
Point* Track::getEnd() {
    return &_end;
}


/**
 * Returns the track's start point
 * @return a pointer to the track's start point
 */
Point* Track::getStart() {
    return &_start;
}


/**
 * Return the track's azimuthal angle (with respect to the x-axis)
 * @return the aximuthal angle
 */
double Track::getPhi() const {
    return _phi;
}


/**
 * Return the track's azimuthal weight
 * @param the track's azimuthal weight
 */
double Track::getWeight() const {
    return _weight;
}


/**
 * Returns a pointer to a segment with a given index or ends program if
 * track does not have the segment
 * @param segment index into the track's segments container
 * @return a pointer to the requested segment
 */
segment* Track::getSegment(int segment) {
	/* Checks to see if segments container contains this segment index */
	if ((int)_segments.size() <= segment)
		return _segments.at(segment);

	/* If track doesn't contain this segment, exits program */
	else
		log_printf(ERROR, "Attempted to retrieve segment s = %d but track only"
				"has %d segments", segment, _segments.size());
	exit(0);
}

/**
 * Return the number of segments along this track
 * @return the number of segments
 */
int Track::getNumSegments() {
	return _segments.size();
}


/**
 * Checks whether a point is contained along this track
 * @param a pointer to the point of interest
 * @return true if the point is on the track, false otherwise
 */
bool Track::contains(Point* point) {

	double m; 		// the slope of the track

	// If the point is outside of the bounds of the start and end points of the track it
	// does not lie on the track
	if (!(((point->getX() <= _start.getX()+1.0E-2 && point->getX() >= _end.getX()-1.0E-2)
		|| (point->getX() >= _start.getX()-1.0E-2 && point->getX() <= _end.getX()+1.0E-2))
		&&
		((point->getY() <= _start.getY()+1.0E-2 && point->getY() >= _end.getY()-1.0E-2)
		|| (point->getY() >= _start.getY()-1.0E-2 && point->getY() <= _end.getY()+1.0E-2))))
		return false;


	// If the track is vertical
	if (fabs(_phi - M_PI / 2) < 1E-10) {
		if (fabs(point->getX() - _start.getX()) < 1E-10)
			return true;
		else
			return false;
	}
	// If the track is not vertical
	else {
		m = sin(_phi) / cos(_phi);

		// Use point-slope formula
		if (fabs(point->getY() - (_start.getY() + m * (point->getX() - _start.getX()))) < 1e-10)
			return true;
		else
			return false;
	}
}


/**
 * Deletes each of this track's segments
 */
void Track::clearSegments() {
	_segments.clear();
}
