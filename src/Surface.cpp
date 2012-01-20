/*
 * Surface.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Surface.h"


/**
 * Default surface constructor
 * @param id the surface id
 * @param type the surface type
 */
Surface::Surface(int id, surfaceType type){
	_id = id;
	_type = type;
}


/**
 * Return the surface's id
 * @return the surface's id
 */
int Surface::getId() const {
    return _id;
}


/**
 * Return the surface's type
 * @return the surface type
 */
surfaceType Surface::getType() const
{
    return _type;
}

/**
 * Destructor
 */
Surface::~Surface() { }


/**
 * Plane constructor
 * @param id the surface id
 * @param A the first coefficient
 * @param B the second coefficient
 * @param C the third coefficient
 */
Plane::Plane(int id, double A, double B, double C): Surface(id, PLANE) {
	_A = A;
	_B = B;
	_C = C;
}

//TODO
std::vector<Surface*> Plane::getNeighborPos() {
    return std::vector<Surface*>();
}


//TODO
std::vector<Surface*> Plane::getNeighborNeg() {
    return std::vector<Surface*>();
}

/**
 * Evaluate a point using the plane's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
double Plane::evaluate(Point* point) {
	double x = point->getX();
	double y = point->getY();
	return (_A * x + _B * y + _C);
}

//TODO
bool Plane::positiveSense(Point* point) {
    return false;
}

/**
 * XPlane constructor for a plane parallel to the x-axis
 * @param id the surface id
 * @param the location of the plane along the y-axis
 */
XPlane::XPlane(int id, double C): Plane(id, 1, 0, -C) {
	_type = XPLANE;
}

//TODO
std::vector<Surface*> XPlane::getNeighborPos() {
    return std::vector<Surface*>();
}


//TODO
std::vector<Surface*> XPlane::getNeighborNeg() {
    return std::vector<Surface*>();
}

//TODO
bool XPlane::positiveSense(Point* point) {
    return false;
}


/**
 * YPlane constructor for a plane parallel to the y-axis
 * @param id the surface id
 * @param the location of the plane along the x-axis
 */
YPlane::YPlane(int id, double C): Plane(id, 1, 0, -C) {
	_type = YPLANE;
}

//TODO
std::vector<Surface*> YPlane::getNeighborPos() {
    return std::vector<Surface*>();
}

//TODO
std::vector<Surface*> YPlane::getNeighborNeg() {
    return std::vector<Surface*>();
}

//TODO
bool YPlane::positiveSense(Point* point) {
    return false;
}

/**
 * Circle constructor
 * @param id the surface id
 * @param x the x-coordinte of the circle center
 * @param y the y-coordinate of the circle center
 * @param radius the radius of the circle
 */
Circle::Circle(int id, double x, double y, double radius): Surface(id, CIRCLE) {
	_A = 1;
	_B = 1;
	_C = -2*x;
	_D = -2*y;
	_E = x*x + y*y - radius*radius;
	_radius = radius;
	center.setX(x);
	center.setY(y);
}

//TODO
std::vector<Surface*> Circle::getNeighborPos() {
    return std::vector<Surface*>();
}


//TODO
std::vector<Surface*> Circle::getNeighborNeg() {
    return std::vector<Surface*>();
}

/**
 * Evaluate a point using the circle's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
double Circle::evaluate(Point* point) {
	double x = point->getX();
	double y = point->getY();
	return (_A * x * x + _B * y * y + _C * x + _D * y + _E);
}


//TODO
bool Circle::positiveSense(Point* point) {
    return false;
}
