/*
 * Surface.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <vector>
#include "Point.h"
#include "Track.h"

#include "log.h"

/* Define for compiler */
class Plane;

/**
 * Surface types
 */
enum surfaceType {
	PLANE,
	CIRCLE,
	XPLANE,
	YPLANE,
	QUADRATIC
};

/**
 * Represents a 2-dimensional quadratics surface
 */
class Surface {
protected:
	int _id;
	surfaceType _type;
public:
	Surface(int id, surfaceType type);
	virtual ~Surface();
	int getId() const;
	surfaceType getType() const;
	virtual std::vector<Surface*> getNeighborPos() =0;
	virtual std::vector<Surface*> getNeighborNeg() =0;
	virtual double evaluate(Point* point) =0;
	virtual int intersection(Track* track, Point* points) =0;
	virtual int intersection(Plane* plane, Point* points) =0;
};

/**
 * Represents a plane in 2D as a Surface subclass
 */
class Plane: public Surface {
private:
	double _A, _B, _C;
	friend class Surface;
	friend class Circle;
public:
	Plane(int id, double A, double B, double C);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
	double evaluate(Point* point);
	int intersection(Track* track, Point* points);
	int intersection(Plane* plane, Point* points);
};

/**
 * Represents a plane parallel to the x-axis as a Plane subclass
 */
class XPlane: public Plane {
private:
	double _A, _B, _C;
public:
	XPlane(int id, double C);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
};

/**
 * Represents a plane parallel to the y-axis as a Plane subclass
 */
class YPlane: public Plane {
private:
	double _A, _B, _C;
public:
	YPlane(int id, double C);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
};

/**
 * Represents a circle as a Surface subclass
 */
class Circle: public Surface {
private:
	Point center;
	double _radius;
	double _A, _B, _C, _D, _E;
	friend class Surface;
	friend class Plane;
public:
	Circle(int id, double x, double y, double radius);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
	double evaluate(Point* point);
	int intersection(Track* track, Point* points);
	int intersection(Plane* plane, Point* points);
};

#endif /* SURFACE_H_ */
