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
#include <sstream>
#include <string>
#include "configurations.h"
#include "log.h"


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
	void setCoords(const double x, const double y);
	double getX() const;
	double getY() const;
	void setX(const double x);
	void setY(const double y);
	double distance(const double x, const double y) const;
	double distance(const Point* point);
	std::string toString();
};

#endif /* POINT_H_ */
