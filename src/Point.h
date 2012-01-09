/*
 * Point.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef POINT_H_
#define POINT_H_

#include <math.h>
#include <iostream>


class Point {
protected:
	float x, y;
public:
	Point();
	Point(const Point& point);
	virtual ~Point();
	void setCoords(float x, float y);
	float getX();
	float getY();
	void setX(float x);
	void setY(float y);
	float distance(float x, float y);
	float distance(Point* point);
};

#endif /* POINT_H_ */
