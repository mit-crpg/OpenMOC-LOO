/*
 * Solver.h
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <utility>
#include "Geometry.h"
#include "Quadrature.h"
#include "Track.h"
#include "TrackGenerator.h"
#include "FlatSourceRegion.h"
#include "configurations.h"
#include "log.h"

class Solver {
private:
	Geometry* _geom;
	Quadrature* _quad;
	FlatSourceRegion* _flat_source_regions;
	Track** _tracks;
	int* _num_tracks;
	int _num_azim;
	int _num_FSRs;
	double _k_eff;
	double _k_eff_old;
public:
	Solver(Geometry* geom, TrackGenerator* track_generator);
	virtual ~Solver();
};

#endif /* SOLVER_H_ */
