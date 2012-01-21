/*
 * Point.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef POINT_H_
#define POINT_H_

#include <math.h>
#include <iostream>
#include "configurations.h"


/**
 * Class to represent a point in 2D
 */
class Point {
private:
	double _x, _y;
public:
	Point();
	Point(const Point& point);
	virtual ~Point();
	void setCoords(double x, double y);
	double getX();
	double getY();
	void setX(double x);
	void setY(double y);
	double distance(double x, double y);
	double distance(Point* point);
};

#endif /* POINT_H_ */
