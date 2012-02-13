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
	_precomputed = false;

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


void Solver::precomputeFactors() {

	Track* curr_track;
	double azim_weight;

	if (_precomputed)
		return;

	/* Precompute the total azimuthal weight for tracks at each polar angle */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];
			azim_weight = curr_track->getAzimuthalWeight();

			for (int p = 0; p < NUM_POLAR_ANGLES; p++) {
				curr_track->setPolarWeight(p, azim_weight*_quad->getMultiple(p));
			}
		}
	}

#ifdef PRECOMPUTE_FACTORS

	segment* curr_seg;

	/* Loop over azimuthal angle, track, segment, polar angle, energy group */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];

			for (int s = 0; s < curr_track->getNumSegments(); s++) {
				curr_seg = curr_track->getSegment(s);

				for (int p = 0; p < NUM_POLAR_ANGLES; p++) {
					for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
						curr_seg->prefactors[p][e] = computePreFactor(curr_seg, e, p);
					}
				}
			}
		}
	}


/* Use hash map */
#else
	log_printf(ERROR, "Table lookup for exponential pre-factors is not yet"
			"implemented. Please set PRECOMPUTE_FACTORS to TRUE in"
			"configurations.h");
#endif

	_precomputed = true;
	return;
}



double Solver::computePreFactor(segment* seg, int energy, int angle) {
	double* sigma_t = seg->_material->getSigmaT();
	double prefactor = 1.0 - exp (-sigma_t[energy] * seg->_length
						/ _quad->getSinTheta(angle));
	return prefactor;
}

