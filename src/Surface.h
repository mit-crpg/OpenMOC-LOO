/*
 * Surface.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include "Point.h"

//using namespace std;

enum surf_type{ CIRCLE, XPLANE, YPLANE, PLANE, QUADRATIC };

class Surface {
protected:
	int uid;
	surf_type type;
	float A, B, C, D, E, F;
public:
	Surface();
	virtual ~Surface();
	float evaluate(Point* point);
};

#endif /* SURFACE_H_ */
