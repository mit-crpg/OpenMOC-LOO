/*
 * TrackGenerator.h
 *
 *  Created on: Jan 23, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include "Point.h"
#include "Track.h"
#include "Geometry.h"

class TrackGenerator {
private:
	int _num_azim;			/* number of azimuthal angles */
	double _spacing;		/* track spacing */
	int* _num_tracks;		/* number of tracks at each angle */
	int* _num_x;			/* number of tracks starting on x-axis */
	int* _num_y;			/* number of tracks starting on y-axis */
	double* _azim_weights;	/* azimuthal weights */
	Track** _tracks;
	Geometry* _geom;
public:
	TrackGenerator(Geometry* geom, int num_azim, double spacing);
	virtual ~TrackGenerator();
    double *getAzimWeights() const;
    int getNumAzim() const;
    int *getNumTracks() const;
    double getSpacing() const;
    Track **getTracks() const;
    void generateTracks();
	void computeEndPoint(Point* start, Point* end,  double phi, double width, double height);
	void makeReflective();
};

#endif /* TRACKGENERATOR_H_ */
