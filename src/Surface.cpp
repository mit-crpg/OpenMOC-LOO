/*
 * Surface.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Surface.h"

/* _n keeps track of the number of surfaces instantiated */
int Surface::_n = 0;


/**
 * Default Surface constructor
 * @param id the surface id
 * @param type the surface type
 * @param boundary this surface's boundary type
 */
Surface::Surface(const int id, const surfaceType type,
					const boundaryType boundary){
	_uid = _n;
	_id = id;
	_type = type;
	_n++;
	_boundary = boundary;
}


/**
 * Surface Destructor
 */
Surface::~Surface() { }


/**
 * Return the surface's uid
 * @return the surface's uid
 */
int Surface::getUid() const {
	return _uid;
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
surfaceType Surface::getType() const {
    return _type;
}


/**
 * Returns the vector of the surface ids on the positive side of this surface
 * @return vector of surface ids
 */
std::vector<Cell*> Surface::getNeighborPos() {
	return _neighbor_pos;
}


/**
 * Returns the vector of the surface ids on the negative side of this surface
 * @return vector of surface id
 */
std::vector<Cell*> Surface::getNeighborNeg() {
	return _neighbor_neg;
}


/**
 * Allocates memory for a certain number of neighbor cells on
 * the positive side of this surface
 * @param size the number of cells
 */
void Surface::setNeighborPosSize(int size) {
	try {
		_neighbor_pos.resize(size);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Could not resize a vector for the positive"
				" neighbor cells for surface id = %d. Backtrace:\n%s",
				_id, e.what());
	}
}


/**
 * Allocates memory for a certain number of neighbor cells on
 * the negative side of this surface
 * @param size the number of cells
 */
void Surface::setNeighborNegSize(int size) {
	try {
		_neighbor_neg.resize(size);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Could not resize a vector for the positive"
				" neighbor cells for surface id = %d. Backtrace:\n%s",
				_id, e.what());
	}
}


/**
 * Sets the id for a neighboring cell on the positive side
 * of this surface
 * @param index the index of the neighbor cell
 * @param cell the cell id
 */
void Surface::setNeighborPos(int index, Cell* cell) {
	_neighbor_pos[index] = cell;
}


/**
 * Sets the id for a neighboring cell on the negative side
 * of this surface
 * @param index the index of the neighbor cell
 * @param cell the cell id
 */
void Surface::setNeighborNeg(int index, Cell* cell) {
	_neighbor_neg[index] = cell;
}


/**
 * Returns the surface's boundary type
 * @return the boundary type (REFLECTIVE or BOUNDARY_NONE)
 */
boundaryType Surface::getBoundary(){
       return _boundary;
}


/**
 * Return true or false if a point is on or off of a surface.
 * @param point pointer to the point of interest
 * @return true (on) or false (off)
 */
bool Surface::onSurface(Point* point) {

	/* Uses a threshold to determine whether the point is on the surface */
	if (abs(evaluate(point)) < ON_SURFACE_THRESH)
		return true;
	else
		return false;
}


/**
 * Return true or false if a localcoord is on or off of a surface
 * @param point pointer to the localcoord of interest
 * @return true (on) or false (off)
 */
bool Surface::onSurface(LocalCoords* coord) {
	return onSurface(coord->getPoint());
}


/**
 * Finds the minimum distance from a point with a given trajectory defined
 * by an angle to this surface. If the trajectory will not intersect the
 * surface, returns INFINITY
 * @param point a pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param intersection a pointer to a point for storing the intersection
 * @return the minimum distance
 */
double Surface::getMinDistance(Point* point, double angle,
									Point* intersection) {

	/* Point array for intersections with this surface */
	Point intersections[2];

	/* Find the intersection point(s) */
	int num_inters = this->intersection(point, angle, intersections);
	double distance = INFINITY;

	/* If there is one intersection point */
	if (num_inters == 1) {
		distance = intersections[0].distance(point);
		intersection->setX(intersections[0].getX());
		intersection->setY(intersections[0].getY());
	}

	/* If there are two intersection points */
	else if (num_inters == 2) {
		double dist1 = intersections[0].distance(point);
		double dist2 = intersections[1].distance(point);

		/* Determine which intersection point is nearest */
		if (dist1 < dist2) {
			distance = dist1;	double getRadius();

			intersection->setX(intersections[0].getX());
			intersection->setY(intersections[0].getY());
		}
		else {
			distance = dist2;
			intersection->setX(intersections[1].getX());
			intersection->setY(intersections[1].getY());
		}
	}

	return distance;
}



/**
 * Plane constructor
 * @param id the surface id
 * @param boundary this surface's boundary type
 * @param A the first coefficient in A * x + B * y + C = 0;
 * @param B the second coefficient in A * x + B * y + C = 0;
 * @param C the third coefficient in A * x + B * y + C = 0;
 */
Plane::Plane(const int id, const boundaryType boundary,
	     const double A, const double B,
	     const double C): Surface(id, PLANE, boundary) {
	_A = A;
	_B = B;
	_C = C;
}


/**
 * Evaluate a point using the plane's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
double Plane::evaluate(const Point* point) const {
	double x = point->getX();
	double y = point->getY();
	return (_A * x + _B * y + _C);
}


/**
 * Converts this Plane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string Plane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = PLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C;

	return string.str();
}



/**
 * Finds the intersection point with this plane from a given point and
 * trajectory defined by an angle
 * @param point pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param points pointer to a point to store the intersection point
 * @return the number of intersection points (0 or 1)
 */
int Plane::intersection(Point* point, double angle, Point* points) {

	double x0 = point->getX();
	double y0 = point->getY();

	int num = 0; 			/* number of intersections */
	double xcurr, ycurr;	/* coordinates of current intersection point */

	/* The track is vertical */
	if ((fabs(angle - (M_PI / 2))) < 1.0e-10) {

		/* The plane is also vertical => no intersections */
		if (_B == 0)
			return 0;

		/* The plane is not vertical */
		else {
			xcurr = x0;
			ycurr = (-_A * x0 - _C) / _B;
			points->setCoords(xcurr, ycurr);
			/* Check that point is in same direction as angle */
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {
		double m = sin(angle) / cos(angle);

		/* The plane and track are parallel, no intersections */
		if (fabs(-_A/_B - m) < 1e-11 && _B != 0)
			return 0;

		else {
			xcurr = -(_B * (y0 - m * x0) + _C)
					/ (_A + _B * m);
			ycurr = y0 + m * (xcurr - x0);
			points->setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;

			return num;
		}
	}
}


/**
 * Finds the intersection points (0 or 1) of a track with a plane.
 * @param track the track of interest
 * @param points an array of two points where intersection points are stored
 * @return the number of intersection points (0 or 1)
 */
int Plane::intersection(Track* track, Point* points) const {

	double x0 = track->getStart()->getX();
	double y0 = track->getStart()->getY();

	int num = 0; 			/* number of intersections */
	double xcurr, ycurr;	/* coordinates of current intersection point */

	/* The track is vertical */
	if ((fabs(track->getPhi() - (M_PI / 2))) < 1.0e-10) {

		/* The plane is also vertical => no intersections */
		if (_B == 0)
			return 0;

		/* The plane is not vertical */
		else {
			xcurr = x0;
			ycurr = (-_A * x0 - _C) / _B;
			points->setCoords(xcurr, ycurr);
			if (track->contains(points))
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {
		double m = sin(track->getPhi()) / cos(track->getPhi());

		/* The plane and track are parallel, no intersections */
		if (fabs(-_A/_B - m) < 1e-5 && _B != 0)
			return 0;

		else {
			xcurr = -(_B * (y0 - m * x0) + _C)
					/ (_A + _B * m);
			ycurr = y0 + m * (xcurr - x0);
			points->setCoords(xcurr, ycurr);
			if (track->contains(points))
				num++;
			return num;
		}
	}
}


/**
 * Finds the intersection points (0 or 1) of a plane with this plane.
 * @param plane the plane of interest
 * @param points an array of two points where intersection points are stored
 * @return the number of intersection points (0 or 1 for a plane)
 */
int Plane::intersection(Plane* plane, Point* points) const {

	double xcurr, ycurr;
	int num = 0;			/* Number of intersection points */

	/* If neither plane is vertical */
	if (_B != 0 && plane->_B != 0) {

		/* If planes are parallel, return 0 */
 		if (_A == plane->_A)
			return num;
 		/* If planes are not parallel, they have an intersection point */
		else {
			xcurr = (plane->_C - _C) / (_A - plane->_A);
			ycurr = -plane->_C - plane->_A * xcurr;
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;
			return num;
		}
	}

	/* If this plane is vertical but parameter plane is non-vertical */
	else if (_B == 0 && plane->_B != 0) {
		xcurr = -_C;
		ycurr = _C * plane->_A - plane->_C;
		points[num].setX(xcurr);
		points[num].setY(ycurr);
		num++;
		return num;
	}

	/* If this plane is non-vertical but parameter plane is vertical */
	else if (_B != 0 && plane->_B == 0) {
		xcurr = -plane->_C;
		ycurr = plane->_C * _A - _C;
		points[num].setX(xcurr);
		points[num].setY(ycurr);
		num++;
		return num;
	}

	/* Both planes are vertical => no intersections points */
	else
		return num;
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double Plane::getXMin(){
	log_printf(ERROR, "Plane::getXMin not implemented");
	return -1.0/0.0;
}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double Plane::getXMax(){
	log_printf(ERROR, "Plane::getXMax not implemented");
	return 1.0/0.0;
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double Plane::getYMin(){
	log_printf(ERROR, "Plane::getYMin not implemented");
	return -1.0/0.0;
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double Plane::getYMax(){
	log_printf(ERROR, "Plane::getYMax not implemented");
	return 1.0/0.0;
}

/**
 * XPlane constructor for a plane parallel to the x-axis
 * @param id the surface id
 * @param boundary this surface's boundary type
 * @param C the location of the plane along the y-axis
 */
XPlane::XPlane(const int id, const boundaryType boundary,
		const double C): Plane(id, boundary, 0, 1, -C) {
	_type = XPLANE;
}


/**
 * Converts this XPlane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string XPlane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = XPLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C;

	return string.str();
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double XPlane::getXMin(){
	return -_C;
}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double XPlane::getXMax(){
	return -_C;
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double XPlane::getYMin(){
	return 1.0/0.0;
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double XPlane::getYMax(){
	return -1.0/0.0;
}


/**
 * YPlane constructor for a plane parallel to the y-axis
 * @param id the surface id
 * @param boundary the surface boundary type
 * @param C the location of the plane along the x-axis
 */
YPlane::YPlane(const int id, const boundaryType boundary,
		const double C): Plane(id, boundary, 1, 0, -C) {
	_type = YPLANE;
}



/**
 * Converts this YPlane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string YPlane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = YPLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C;

	return string.str();
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double YPlane::getXMin(){
	return 1.0/0.0;
}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double YPlane::getXMax(){
	return -1.0/0.0;
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double YPlane::getYMin(){
	return -_C;
}

/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double YPlane::getYMax(){
	return -_C;
}


/**
 * Circle constructor
 * @param id the surface id
 * @param boundary the surface boundary type
 * @param x the x-coordinte of the circle center
 * @param y the y-coordinate of the circle center
 * @param radius the radius of the circle
 */
Circle::Circle(const int id, const boundaryType boundary, const double x,
		const double y, const double radius): Surface(id, CIRCLE, boundary) {
	_A = 1;
	_B = 1;
	_C = -2*x;
	_D = -2*y;
	_E = x*x + y*y - radius*radius;
	_radius = radius;
	center.setX(x);
	center.setY(y);
}

/**
 * Return the radius of the circle
 * @return the radius of the circle
 */
double Circle::getRadius() {
	return this->_radius;
}



/**
 * Evaluate a point using the circle's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
double Circle::evaluate(const Point* point) const {
	double x = point->getX();
	double y = point->getY();
	return (_A * x * x + _B * y * y + _C * x + _D * y + _E);
}



/**
 * Finds the intersection point with this circle from a given point and
 * trajectory defined by an angle (0, 1, or 2 points)
 * @param point pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param points pointer to a an array of points to store intersection points
 * @return the number of intersection points (0 or 1)
 */
int Circle::intersection(Point* point, double angle, Point* points) {

	double x0 = point->getX();
	double y0 = point->getY();
	double xcurr, ycurr;
	int num = 0;			/* Number of intersection points */
	double a, b, c, q, discr;

	/* If the track is vertical */
	if ((fabs(angle - (M_PI / 2))) < 1.0e-10) {
		/* Solve for where the line x = x0 and the surface F(x,y) intersect
		 * Find the y where F(x0, y) = 0
		 * Substitute x0 into F(x,y) and rearrange to put in
		 * the form of the quadratic formula: ay^2 + by + c = 0
		 */
		a = _B * _B;
		b = _D;
		c = _A * x0 * x0 + _C * x0 + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = x0;
			ycurr = -b / (2*a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = x0;
			ycurr = (-b + sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;

			xcurr = x0;
			ycurr = (-b - sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {
		/*Solve for where the line y-y0 = m*(x-x0) and the surface F(x,y) intersect
		 * Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
		 * Substitute the point-slope formula for y into F(x,y) and rearrange to put in
		 * the form of the quadratic formula: ax^2 + bx + c = 0
		 * double m = sin(track->getPhi()) / cos(track->getPhi());
		 */
		double m = sin(angle) / cos(angle);
		q = y0 - m * x0;
		a = _A + _B * _B * m * m;
		b = 2 * _B * m * q + _C + _D * m;
		c = _B * q * q + _D * q + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = -b / (2*a);
			ycurr = y0 + m * (points[0].getX() - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = (-b + sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0) {
				num++;
			}
			else if (angle > M_PI && ycurr < y0) {
				num++;
			}

			xcurr = (-b - sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0) {
				num++;
			}
			else if (angle > M_PI && ycurr < y0) {
				num++;
			}

			return num;
		}
	}
}


/**
 * Finds the intersection points (0, 1, or 2) of a track with a circle
 * @param track the track of interest
 * @param points an array of two points where the intersection points are stored
 * @return the number of intersection points (0, 1 or 2 for a circle)
 */
int Circle::intersection(Track* track, Point* points) const {

	double x0 = track->getStart()->getX();
	double y0 = track->getStart()->getY();
	double xcurr, ycurr;
	int num = 0;			/* Number of intersection points */	double getRadius();

	double a, b, c, q, discr;

	/* If the track is vertical */
	if ((fabs(track->getPhi() - (M_PI / 2))) < 1.0e-10) {
		/* Solve for where the line x = x0 and the surface F(x,y) intersect
		 * Find the y where F(x0, y) = 0			if (track->contains(points))
				num++;
		 *
		 * Substitute x0 into F(x,y) and rearrange to put in
		 * the form of the quadratic formula: ay^2 + by + c = 0
		 */
		a = _B * _B;
		b = _D;
		c = _A * x0 * x0 + _C * x0 + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = x0;
			ycurr = -b / (2*a);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;

			return num;
		}

		/* There are two intersections */
		else {
			xcurr = x0;
			ycurr = (-b + sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;

			xcurr = x0;
			ycurr = (-b - sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {

		/*Solve for where the line y-y0 = m*(x-x0) and the surface F(x,y)
		 * intersect
		 * Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
		 * Substitute the point-slope formula for y into F(x,y) and rearrange
		 * to put in the form of the quadratic formula: ax^2 + bx + c = 0
		 * double m = sin(track->getPhi()) / cos(track->getPhi());
		 */
		double m = sin(track->getPhi()) / cos(track->getPhi());
		q = y0 - m * x0;
		a = _A + _B * _B * m * m;
		b = 2 * _B * m * q + _C + _D * m;
		c = _B * q * q + _D * q + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = -b / (2*a);
			ycurr = y0 + m * (points[0].getX() - x0);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;
			return num;
		}


		/* There are two intersections */
		else {
			xcurr = (-b + sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;

			xcurr = (-b - sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (track->contains(&points[num]))
				num++;

			return num;
		}
	}
}


/**
 * Finds the intersection points (0, 1, or 2) of a plane with a circle
 * @param plane the plane of interest
 * @param points an array of two points where intersection points are stored
 * @return the number of intersection points (0, 1 or 2 for a circle)
 */
int Circle::intersection(Plane* plane, Point* points) const {

	float xcurr, ycurr;
	int num = 0;
	float a, b, c, discr;

	/* If the plane is vertical */
	if (plane->_B == 0) {

		/* Solve for where the plane F(x,y) = x + F1 and
		 * the circle F(x,y) = x^2 + y^2 + D2*x + E2*y + F2 intersect
		 * Set F(x,y) for the plane equal to zero, rearrange for x,
		 * and substitute into the F(x,y) for the circle.
		 * Rearrange to put in the form of the quadratic formula:
		 * ay^2 + by + c = 0
		 * Solve for y and find x from the equation for the plane
		 */
		a = 1;
		b = _D;
		c = _E - plane->_C * this->_C + plane->_C * plane->_C;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = -plane->_C;
			ycurr = -b / (2*a);
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = -plane->_C;
			ycurr = (-b + sqrt(discr)) / (2*a);
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;

			ycurr = (-b - sqrt(discr)) / (2*a);
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;

			return num;
		}
	}

	/* If the plane isn't vertical */
	else {
		/* Solve for where the plane F(x,y) = D1*x + y + F1 and
		 * the circle F(x,y) = x^2 + y^2 + D2*x + E2*y + F2 intersect
		 * Set F(x,y) for the plane equal to zero, rearrange for y,
		 * and substitute into the F(x,y) for the circle.
		 * Rearrange to put in the form of the quadratic formula:
		 * ax^2 + bx + c = 0	double getRadius();
		 *
		 * Solve for x and find y from the equation for the plane
		 */
		a = 1 + plane->_A * plane->_A;
		b = 2 * plane->_C * plane->_A + this->_C - this->_D * plane->_A;
		c = plane->_C * plane->_C + this->_E - this->_D * plane->_C;

		discr = b * b - 4 * a * c;

		/* There are no intersections */
		if (discr < 0) {
			return 0;
		}

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = -b / (2*a);
			ycurr = -plane->_C - plane->_A * xcurr;
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = (-b + sqrt(discr)) / (2*a);
			ycurr = -plane->_C - plane->_A * xcurr;
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;

			xcurr = (-b - sqrt(discr)) / (2*a);
			ycurr = -plane->_C - plane->_A * xcurr;
			points[num].setX(xcurr);
			points[num].setY(ycurr);
			num++;
			return num;
		}
	}
}


/**
 * Converts this Plane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string Circle::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = CIRCLE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D
			<< ", E = " << _E;

	return string.str();
}

double Circle::getXMin(){
	log_printf(ERROR, "Circle::getXMin not implemented");
	return -1.0/0.0;
}

double Circle::getXMax(){
	log_printf(ERROR, "Circle::getXMax not implemented");
	return 1.0/0.0;
}

double Circle::getYMin(){
	log_printf(ERROR, "Circle::getYMin not implemented");
	return -1.0/0.0;
}

double Circle::getYMax(){
	log_printf(ERROR, "Circle::getYMax not implemented");
	return 1.0/0.0;
}


