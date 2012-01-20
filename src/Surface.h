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

/**
 * Surface types
 */
enum surfaceType {
	CIRCLE, XPLANE, YPLANE, PLANE, QUADRATIC
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
	virtual std::vector<Surface*> getNeighborPos() =0;
	virtual std::vector<Surface*> getNeighborNeg() =0;
	virtual double evaluate(Point* point) =0;
	virtual bool positiveSense(Point* point) =0;
    int getId() const;
    surfaceType getType() const;
//	virtual intersection(Track* track) =0;
//	virtual intersection(Plane* plane) =0;
};

/**
 * Represents a plane in 2D as a Surface subclass
 */
class Plane: public Surface {
private:
	double _A, _B, _C;
	friend class Surface;
public:
	Plane(int id, double A, double B, double C);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
	double evaluate(Point* point);
	bool positiveSense(Point* point);
};

/**
 * Represents a plane parallel to the x-axis as a Plane subclass
 */
class XPlane: public Plane {
private:
	double _A, _B, _C;
//	friend class Plane;
//	friend class Surface;
public:
	XPlane(int id, double C);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
	bool positiveSense(Point* point);
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
	bool positiveSense(Point* point);
//	friend class Plane;
//	friend class Surface;
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
public:
	Circle(int id, double x, double y, double radius);
	std::vector<Surface*> getNeighborPos();
	std::vector<Surface*> getNeighborNeg();
	double evaluate(Point* point);
	bool positiveSense(Point* point);
};

#endif /* SURFACE_H_ */
