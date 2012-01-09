/*
 * Point.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
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
 * Initialize a point
 */
void Point::setCoords(float x, float y) {
	this->x = x;
	this->y = y;
}


/**
 * Get the x-coordinate
 */
float Point::getX() {
	return this->x;
}


/**
 * Get the y-coordinate
 */
float Point::getY() {
	return this->y;
}


/**
 * Set the x-coordinate
 */
void Point::setX(float x) {
	this->x = x;
}


/**
 * Set the y-coordinate
 */
void Point::setY(float y) {
	this->y = y;
}


/**
 * Compute the distance from this point to a point of interest
 */
float Point::distance(float x, float y) {
	float deltax = this->x - x;
	float deltay = this->y - y;
	return sqrt(deltax*deltax + deltay*deltay);
}



/**
 * Compute the distance from this point to a point of interest
 */
float Point::distance(Point* point) {
	float deltax = this->x - point->x;
	float deltay = this->y - point->y;
	return sqrt(deltax*deltax + deltay*deltay);
}
