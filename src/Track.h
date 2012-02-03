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
#include <string>
#include "Point.h"
#include "log.h"

/* Represent a segment along a given track */
struct segment {
	double _length;
	int _region_id;
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
	void setValues(const double start_x, const double start_y,
			const double end_x, const double end_y, const double phi);
    void setWeight(const double weight);
    void setPhi(const double phi);
    void setReflIn(const bool refl_in);
    void setReflOut(const bool refl_out);
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
	std::string toString();
};

#endif /* TRACK_H_ */
