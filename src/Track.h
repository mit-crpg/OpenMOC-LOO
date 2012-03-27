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
#include "configurations.h"

#if USE_OPENMP
	#include <omp.h>
#endif

/* Represent a segment along a given track */
struct segment {
	double _length;
	Material* _material;
	int _region_id;
#if STORE_PREFACTORS
	double _prefactors[NUM_ENERGY_GROUPS][NUM_POLAR_ANGLES];
#endif
};

class Track {
private:
	Point _start;
	Point _end;
	double _phi;
	double _azim_weight;
	double _polar_weights[NUM_POLAR_ANGLES];
	double _polar_fluxes[2 * GRP_TIMES_ANG];
	std::vector<segment*> _segments;
	Track *_track_in, *_track_out;
	bool _refl_in, _refl_out;
#if USE_OPENMP
	omp_lock_t _flux_lock;
#endif
public:
	Track();
	virtual ~Track();
	void setValues(const double start_x, const double start_y,
			const double end_x, const double end_y, const double phi);
    void setAzimuthalWeight(const double azim_weight);
    void setPolarWeight(const int angle, double polar_weight);
    void setPolarFluxes(bool direction, int start_index, double* polar_fluxes);
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
    double* getPolarFluxes();
	segment* getSegment(int s);
	std::vector<segment*> getSegments();
	int getNumSegments();
    Track *getTrackIn() const;
    Track *getTrackOut() const;
    bool isReflIn() const;
    bool isReflOut() const;

    void normalizeFluxes(double factor);
    bool contains(Point* point);
	void addSegment(segment* segment);
	void clearSegments();
	std::string toString();
};

#endif /* TRACK_H_ */
