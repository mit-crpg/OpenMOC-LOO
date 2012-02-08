/*
 * Solver.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Solver.h"

Solver::Solver(Geometry* geom, TrackGenerator* track_generator) {
	_geom = geom;
	_num_FSRs = geom->getNumFSRs();
	_tracks = track_generator->getTracks();
	_num_tracks = track_generator->getNumTracks();
	_num_azim = track_generator->getNumAzim();

	try{
		_flat_source_regions = new FlatSourceRegion[_num_FSRs];
	}
	catch(std::exception &e) {
		log_printf(ERROR, "Could not allocate memory for the solver's flat "
					"source region array. Backtrace:%s", e.what());
	}

}

Solver::~Solver() {
	delete _flat_source_regions;
}




