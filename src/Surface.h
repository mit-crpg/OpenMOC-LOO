/*
 * Surface.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <vector>
#include <sstream>
#include <string>
#include "Point.h"
#include "Track.h"
#include "Cell.h"

#include "log.h"

/* Define for compiler */
class Plane;
class Cell;

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

enum boundaryType {
	BOUNDARY_NONE,
	REFLECTIVE
};



/**
 * Represents a 2-dimensional quadratics surface
 */
class Surface {
protected:
	static int _n;				/* Counts the number of surfaces */
	int _uid;				/* monotonically increasing id based on n */
	int _id;
	surfaceType _type;
	boundaryType _boundary;
	std::vector<Cell*> _neighbor_pos;
	std::vector<Cell*> _neighbor_neg;

public:
	Surface(const int id, const surfaceType type, const boundaryType boundary);
	virtual ~Surface();
	int getUid() const;
	int getId() const;
	surfaceType getType() const;
//	std::vector<int> getNeighborPos();
//	std::vector<int> getNeighborNeg();
	std::vector<Cell*> getNeighborPos();
	std::vector<Cell*> getNeighborNeg();

	void setNeighborPosSize(int size);
	void setNeighborNegSize(int size);
	void setNeighborPos(int index, Cell* cell);
	void setNeighborNeg(int index, Cell* cell);
	boundaryType getBoundary();
	virtual double evaluate(const Point* point) const =0;
	virtual int intersection(Track* track, Point* points) const =0;
	virtual int intersection(Plane* plane, Point* points) const =0;
	virtual std::string toString() =0;
	virtual double getXMin() =0;
	virtual double getXMax() =0;
	virtual double getYMin() =0;
	virtual double getYMax() =0;
};

/**
 * Represents a plane in 2D as a Surface subclass
 */
class Plane: public Surface {
protected:
	double _A, _B, _C;
	friend class Surface;
	friend class Circle;
public:
	Plane(const int id, const boundaryType boundary, const double A, const double B, const double C);
	double evaluate(const Point* point) const;
	int intersection(Track* track, Point* points) const;
	int intersection(Plane* plane, Point* points) const;
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();
};

/**
 * Represents a plane parallel to the x-axis as a Plane subclass
 */
class XPlane: public Plane {
private:
public:
	XPlane(const int id, const boundaryType boundary, const double C);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();
};

/**
 * Represents a plane parallel to the y-axis as a Plane subclass
 */
class YPlane: public Plane {
private:
public:
	YPlane(const int id, const boundaryType boundary, const double C);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();
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
	Circle(const int id, const boundaryType boundary, const double x,
				const double y, const double radius);
	double evaluate(const Point* point) const;
	int intersection(Track* track, Point* points) const;
	int intersection(Plane* plane, Point* points) const;
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();
	
};

#endif /* SURFACE_H_ */
