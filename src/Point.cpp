/*
 * Point.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Point.h"


/**
 * Default constructor
 */
Point::Point() { }


/**
 * Destructor
 */
Point::~Point() { }


/**
 * Initializes a point
 * @param x x-coordinate
 * @param y y-coordinate
 */
void Point::setCoords(const double x, const double y) {
	_x = x;
	_y = y;
}


/**
 * Returns this point's x-coordinate
 * @return the x-coordinate
 */
double Point::getX() const {
	return _x;
}


/**
 * Returns this point's y-coordinate
 * @return the y-coordinate
 */
double Point::getY() const {
	return _y;
}


/**
 * Set the point's x-coordinate
 * @param x the new x-coordinate
 */
void Point::setX(const double x) {
	_x = x;
}


/**
 * Set the point's y-coordinate
 * @param y the new y-coordinate
 */
void Point::setY(const double y) {
	_y = y;
}


/**
 * Compute the distance from this point to a point of interest
 * @param x the x-coordinate of the point of interest
 * @param y the y-coordinate of the point of interest
 * @return distance to the point of interest
 */
double Point::distance(const double x, const double y) const {
	double deltax = _x - x;
	double deltay = _y - y;
	return sqrt(deltax*deltax + deltay*deltay);
}



/**
 * Compute the distance from this point to a point of interest
 * @param point the point of interest
 * @return distance to the point of interest
 */
double Point::distance(const Point* point) {
	double deltax = _x - point->_x;
	double deltay = _y - point->_y;
	return sqrt(deltax*deltax + deltay*deltay);
}
