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
#include "Material.h"
#include "log.h"

/* Represent a segment along a given track */
struct segment {
	double _length;
	Material* _material;
	int _region_id;
#ifdef PRECOMPUTE_FACTORS
	double prefactors[NUM_POLAR_ANGLES][NUM_ENERGY_GROUPS];
#endif
};


class Track {
private:
	Point _start;
	Point _end;
	double _phi;
	double _azim_weight;
	double _polar_weights[NUM_POLAR_ANGLES];
	std::vector<segment*> _segments;
	Track *_track_in, *_track_out;
	bool _refl_in, _refl_out;
public:
	Track();
	virtual ~Track();
	void setValues(const double start_x, const double start_y,
			const double end_x, const double end_y, const double phi);
    void setAzimuthalWeight(const double azim_weight);
    void setPolarWeight(const int angle, double polar_weight);
    void setPhi(const double phi);
    void setReflIn(const bool refl_in);
    void setReflOut(const bool refl_out);
    void setTrackIn(Track *track_in);
    void setTrackOut(Track *track_out);
    Point* getEnd();
    Point* getStart();
    double getPhi() const;
    double getAzimuthalWeight() const;
    double* getPolarWeights();
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
