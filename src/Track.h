/*
 * Track.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef TRACK_H_
#define TRACK_H_

#include <vector>
#include "Point.h"

class Track {
private:
	Point _start;
	Point _end;
	double _phi;
	double _weight;
	std::vector<segment> _segments;
public:
	Track(double start_x, double start_y, double end_x,
			double end_y, double phi);
	virtual ~Track();
    Point getEnd() const;
    double getPhi() const;
    Point getStart() const;
    double getWeight() const;
};



struct segment {
	double length;
	int region;
};

#endif /* TRACK_H_ */
