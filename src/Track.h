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
#include <stdlib.h>
#include "Point.h"

struct segment {
	double length;
	int region;
};


class Track {
private:
	Point _start;
	Point _end;
	double _phi;
	double _weight;
	std::vector<segment*> _segments;
public:
	Track(double start_x, double start_y, double end_x,
			double end_y, double phi);
	virtual ~Track();
    void setWeight(double weight);
    Point* getEnd();
    Point* getStart();
    double getPhi() const;
    double getWeight() const;
	segment* getSegment(int s);
	int getNumSegments();
	bool contains(Point* point);
	void addSegment(segment* segment);
	void clearSegments();
};



#endif /* TRACK_H_ */
