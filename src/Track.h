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
#include "log.h"

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
	Track *_track_in, *_track_out;
	bool _refl_in, _refl_out;
public:
	Track();
	virtual ~Track();
	void setValues(double start_x, double start_y, double end_x,
			double end_y, double phi);
    void setWeight(double weight);
    void setReflIn(bool refl_in);
    void setReflOut(bool refl_out);
    void setTrackIn(Track *track_in);
    void setTrackOut(Track *track_out);
    Point* getEnd();
    Point* getStart();
    double getPhi() const;
    double getWeight() const;
	segment* getSegment(int s);
	int getNumSegments();
    Track *getTrackIn() const;
    Track *getTrackOut() const;
    bool isReflIn() const;
    bool isReflOut() const;
	bool contains(Point* point);
	void addSegment(segment* segment);
	void clearSegments();

};



#endif /* TRACK_H_ */
