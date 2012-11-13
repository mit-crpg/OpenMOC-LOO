/*
 * Solver.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Solver.h"


/**
 * Solver constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the trackgenerator
 */
Solver::Solver(Geometry* geom, TrackGenerator* track_generator,
								Plotter* plotter, Cmfd* cmfd, bool updateFlux) {
	_geom = geom;
	_quad = new Quadrature(TABUCHI);
	_num_FSRs = geom->getNumFSRs();
	_tracks = track_generator->getTracks();
	_num_tracks = track_generator->getNumTracks();
	_num_azim = track_generator->getNumAzim();
	_plotter = plotter;
	_update_flux = updateFlux;
	_cmfd = cmfd;
	try{
		_flat_source_regions = new FlatSourceRegion[_num_FSRs];
		_FSRs_to_powers = new double[_num_FSRs];
		_FSRs_to_pin_powers = new double[_num_FSRs];

		for (int e = 0; e <= NUM_ENERGY_GROUPS; e++) {
			_FSRs_to_fluxes[e] = new double[_num_FSRs];
			_FSRs_to_absorption[e] = new double[_num_FSRs];
			_FSRs_to_pin_absorption[e] = new double[_num_FSRs];
		}
	}
	catch(std::exception &e) {
		log_printf(ERROR, "Could not allocate memory for the solver's flat "
					"source region array. Backtrace:%s", e.what());
	}

	/* Pre-compute exponential pre-factors */
	precomputeFactors();
	initializeFSRs();
}


/**
 * Solver destructor deletes flat source regions array
 */
Solver::~Solver() {

	delete [] _flat_source_regions;
	delete [] _FSRs_to_powers;
	delete [] _FSRs_to_pin_powers;
	delete _quad;

	for (int e = 0; e <= NUM_ENERGY_GROUPS; e++)
		delete [] _FSRs_to_fluxes[e];

#if !STORE_PREFACTORS
	delete [] _pre_factor_array;
#endif

}



/**
 * Pre-computes exponential pre-factors for each segment of each track for
 * each polar angle. This method will store each pre-factor in an array inside
 * each segment if STORE_PREFACTORS is set to true inside the configurations.h
 * file. If it is not set to true then a hashmap will be generated which will
 * contain values of the pre-factor at for specific segment lengths (the keys
 * into the hashmap).
 */
void Solver::precomputeFactors() {

	log_printf(INFO, "Pre-computing exponential pre-factors...");

	Track* curr_track;
	double azim_weight;

	/* Precompute the total azimuthal weight for tracks at each polar angle */
	#if USE_OPENMP
	#pragma omp parallel for private(curr_track, azim_weight)
	#endif
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];
			azim_weight = curr_track->getAzimuthalWeight();

			for (int p = 0; p < NUM_POLAR_ANGLES; p++)
				curr_track->setPolarWeight(p, azim_weight*_quad->getMultiple(p) * FOUR_PI);
		}
	}


/*Store pre-factors inside each segment */
#if STORE_PREFACTORS

	log_printf(INFO, "Pre-factors will be stored inside each segment...");

	segment* curr_seg;

	/* Loop over azimuthal angle, track, segment, polar angle, energy group */
	#if USE_OPENMP
	#pragma omp parallel for private(curr_track, curr_seg)
	#endif
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];

			for (int s = 0; s < curr_track->getNumSegments(); s++) {
				curr_seg = curr_track->getSegment(s);

				for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
					for (int p = 0; p < NUM_POLAR_ANGLES; p++) {
						curr_seg->_prefactors[e][p] = computePreFactor(curr_seg, e, p);
					}
				}
			}
		}
	}



/* Use hash map */
#else

	/* make pre factor array based on table look up with linear interpolation */

	log_printf(NORMAL, "Making Prefactor array...");

	/* set size of prefactor array */
	int num_array_values = 10 * sqrt(1 / (8 * KEFF_CONVERG_THRESH));
	_pre_factor_spacing = 10.0 / num_array_values;
	_pre_factor_array_size = 2 * NUM_POLAR_ANGLES * num_array_values;
	_pre_factor_max_index = _pre_factor_array_size - 2*NUM_POLAR_ANGLES - 1;

	log_printf(DEBUG, "prefactor array size: %i, max index: %i", _pre_factor_array_size, _pre_factor_max_index);

	/* allocate arrays */
	_pre_factor_array = new double[_pre_factor_array_size];

	double expon;
	double intercept;
	double slope;

	/* Create prefactor array */
	for (int i = 0; i < num_array_values; i ++){
		for (int j = 0; j < NUM_POLAR_ANGLES; j++){
			expon = exp(- (i * _pre_factor_spacing) / _quad->getSinTheta(j));
			slope = - expon / _quad->getSinTheta(j);
			intercept = expon * (1 + (i * _pre_factor_spacing) / _quad->getSinTheta(j));
			_pre_factor_array[2 * NUM_POLAR_ANGLES * i + 2 * j] = slope;
			_pre_factor_array[2 * NUM_POLAR_ANGLES * i + 2 * j + 1] = intercept;
		}
	}

#endif

	return;
}


/**
 * Function to compute the exponential prefactor for the transport equation for
 * a given segment
 * @param seg pointer to a segment
 * @param energy energy group index
 * @param angle polar angle index
 * @return the pre-factor
 */
double Solver::computePreFactor(segment* seg, int energy, int angle) {
	double* sigma_t = seg->_material->getSigmaT();
	double prefactor = 1.0 - exp (-sigma_t[energy] * seg->_length
						/ _quad->getSinTheta(angle));
	return prefactor;
}


/**
 * Compute the ratio of source / sigma_t for each energy group in each flat
 * source region for efficient fixed source iteration
 */
void Solver::computeRatios() {

	#if USE_OPENMP
	#pragma omp parallel for
	#endif
	for (int i = 0; i < _num_FSRs; i++)
		_flat_source_regions[i].computeRatios();

	return;
}


/**
 * Initializes each of the FlatSourceRegion objects inside the solver's
 * array of FSRs. This includes assigning each one a unique, monotonically
 * increasing id, setting the material for each FSR, and assigning a volume
 * based on the cumulative length of all of the segments inside the FSR.
 */
void Solver::initializeFSRs() {

	log_printf(NORMAL, "Initializing FSRs...");

	CellBasic* cell;
	Material* material;
	Universe* univ_zero = _geom->getUniverse(0);
	Track* track;
	segment* seg;
	FlatSourceRegion* fsr;

	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Loop over azimuthal angle, track and segment */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			track = &_tracks[i][j];

			for (int s = 0; s < track->getNumSegments(); s++) {
				seg = track->getSegment(s);
				fsr =&_flat_source_regions[seg->_region_id];
				fsr->incrementVolume(seg->_length * track->getAzimuthalWeight());
			}
		}
	}

	/* Loop over all FSRs */
	#if USE_OPENMP
	#pragma omp parallel for private(cell, material)
	#endif
	for (int r = 0; r < _num_FSRs; r++) {
		/* Set the id */
		_flat_source_regions[r].setId(r);

		/* Get the cell corresponding to this FSR from the geometry */
		cell = static_cast<CellBasic*>(_geom->findCell(univ_zero, r));

		/* Get the cell's material and assign it to the FSR */
		material = _geom->getMaterial(cell->getMaterial());
		_flat_source_regions[r].setMaterial(material);

		log_printf(INFO, "FSR id = %d has cell id = %d and material id = %d "
				"and volume = %f", r, cell->getId(), material->getId(),
				_flat_source_regions[r].getVolume());
	}

	return;
}


/**
 * Zero each track's incoming and outgoing polar fluxes
 */
void Solver::zeroTrackFluxes() {

	log_printf(INFO, "Setting all track polar fluxes to zero...");

	double* polar_fluxes;

	/* Loop over azimuthal angle, track, polar angle, energy group
	 * and set each track's incoming and outgoing flux to zero */
	#if USE_OPENMP
	#pragma omp parallel for private(polar_fluxes)
	#endif
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			polar_fluxes = _tracks[i][j].getPolarFluxes();

			for (int i = 0; i < GRP_TIMES_ANG * 2; i++)
				polar_fluxes[i] = 0.0;
		}
	}
}


/**
 * Set the scalar flux for each energy group inside each FSR
 * to unity
 */
void Solver::oneFSRFluxes() {

	log_printf(INFO, "Setting all FSR scalar fluxes to unity...");
	FlatSourceRegion* fsr;

	/* Loop over all FSRs and energy groups */
	#if USE_OPENMP
	#pragma omp parallel for private(fsr)
	#endif
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_flat_source_regions[r];
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			fsr->setFlux(e, 1.0);
	}

	return;
}


/**
 * Set the scalar flux for each energy group inside each FSR to zero
 */
void Solver::zeroFSRFluxes() {

	log_printf(INFO, "Setting all FSR scalar fluxes to zero...");
	FlatSourceRegion* fsr;

	/* Loop over all FSRs and energy groups */
	#if USE_OPENMP
	#pragma omp for private(fsr)
	#endif
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_flat_source_regions[r];
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			fsr->setFlux(e, 0.0);
	}

	return;
}

/**
 * Set all MeshCell currents, crossings, and weights to zero.
 */
void Solver::zeroMeshCells() {

	/* get mesh */
	Mesh* mesh = _geom->getMesh();
	MeshCell* meshCell;

	/* loop over mesh cells */
	for (int i = 0; i < mesh->getCellHeight()*mesh->getCellWidth(); i++){
		meshCell = mesh->getCells(i);

		/* loop over mesh surfaces in mesh cell */
		for (int surface = 0; surface < 8; surface++){

			/* set number of crossings to 0 */
			meshCell->getMeshSurfaces(surface)->setCrossings(0);

			/* loop over azimuthal angles */
			for (int azim = 0; azim < _num_azim; azim++){

				/* set surface weight to zero */
				meshCell->getMeshSurfaces(surface)->setWeight(0, azim);

				/* loop over energy groups */
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){

					/* set mesh cell fluxes to 0 */
					meshCell->setOldFlux(0, group);
					meshCell->setNewFlux(0, group);

					/* set current to zero */
					meshCell->getMeshSurfaces(surface)->setCurrent(0, group, azim);
				}
			}
		}
	}
}


/**
 * Compute k_eff from the new and old source and the value of k_eff from
 * the previous iteration
 */
void Solver::updateKeff() {

	double tot_abs = 0.0;
	double tot_fission = 0.0;
	double abs = 0;
	double fission;
	double* sigma_a;
	double* nu_sigma_f;
	double* flux;
	Material* material;
	FlatSourceRegion* fsr;

#if USE_OPENMP
#pragma omp parallel shared(tot_abs, tot_fission)
{
	#pragma omp for private(fsr, material, sigma_a, nu_sigma_f, flux, abs, fission)
	#endif
	for (int r = 0; r < _num_FSRs; r++) {
		abs = 0;
		fission = 0;
		fsr = &_flat_source_regions[r];
		material = fsr->getMaterial();
		sigma_a = material->getSigmaA();
		nu_sigma_f = material->getNuSigmaF();
		flux = fsr->getFlux();

		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			abs += sigma_a[e] * flux[e] * fsr->getVolume();
			fission += nu_sigma_f[e] * flux[e] * fsr->getVolume();
		}

			#if USE_OPENMP
			#pragma omp atomic
			#endif
			tot_abs += abs;

			#if USE_OPENMP
			#pragma omp atomic
			#endif
			tot_fission += fission;
	}

#if USE_OPENMP
}
#endif

	_k_eff = tot_fission/tot_abs;
	log_printf(INFO, "Computed k_eff = %f", _k_eff);

	return;
}


/**
 * Return an array indexed by FSR ids which contains the corresponding
 * fluxes for each FSR
 * @return an array map of FSR to fluxes
 */
double** Solver::getFSRtoFluxMap() {
	return _FSRs_to_fluxes;
}


/**
 * Checks that each flat source region has at least one segment within it
 * and if not, exits the program with an error message
 */
void Solver::checkTrackSpacing() {

	int* FSR_segment_tallies = new int[_num_FSRs];
	Track* track;
	std::vector<segment*> segments;
	segment* segment;
	int num_segments;
	Cell* cell;

	/* Set each tally to zero to begin with */
	for (int i=0; i < _num_FSRs; i++)
		FSR_segment_tallies[i] = 0;


	/* Iterate over all azimuthal angles, all tracks, and all segments
	 * and tally each segment in the corresponding FSR */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			track = &_tracks[i][j];
			segments = track->getSegments();
			num_segments = track->getNumSegments();

			for (int s = 0; s < num_segments; s++) {
				segment = segments.at(s);
				FSR_segment_tallies[segment->_region_id]++;
			}
		}
	}


	/* Loop over all FSRs and if one FSR does not have tracks in it, print
	 * error message to the screen and exit program */
	for (int i=0; i < _num_FSRs; i++) {
		if (FSR_segment_tallies[i] == 0) {
			cell = _geom->findCell(i);
			log_printf(ERROR, "No tracks were tallied inside FSR id = %d which "
					"is cell id = %d. Please reduce your track spacing,"
					" increase the number of azimuthal angles, or increase the"
					" size of the flat source regions", i, cell->getId());
		}
	}

	delete [] FSR_segment_tallies;
}


/**
 * Compute the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates
 */
void Solver::computePinPowers() {

	log_printf(NORMAL, "Computing pin powers...");

	FlatSourceRegion* fsr;
	double tot_pin_power = 0;
	double avg_pin_power = 0;
	double num_nonzero_pins = 0;
	double curr_pin_power = 0;
	double prev_pin_power = 0;

	/* create BitMaps for plotting */
	BitMap<int>* bitMapFSR = new BitMap<int>;
	BitMap<float>* bitMap = new BitMap<float>;
	bitMapFSR->pixel_x = _plotter->getBitLengthX();
	bitMapFSR->pixel_y = _plotter->getBitLengthY();
	bitMap->pixel_x = _plotter->getBitLengthX();
	bitMap->pixel_y = _plotter->getBitLengthY();
	initialize(bitMapFSR);
	initialize(bitMap);
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMapFSR->color_type = RANDOM;
	bitMap->color_type = SCALED;

	/* make FSR BitMap */
	_plotter->makeFSRMap(bitMapFSR->pixels);

	/* Loop over all FSRs and compute the fission rate*/
	for (int i=0; i < _num_FSRs; i++) {
		fsr = &_flat_source_regions[i];
		_FSRs_to_powers[i] = fsr->computeFissionRate();
	}

	/* Compute the pin powers by adding up the powers of FSRs in each
	 * lattice cell, saving lattice cell powers to files, and saving the
	 * pin power corresponding to each FSR id in FSR_to_pin_powers */
	_geom->computePinPowers(_FSRs_to_powers, _FSRs_to_pin_powers);


	/* Compute the total power based by accumulating the power of each unique
	 * pin with a nonzero power */
	for (int i=0; i < _num_FSRs; i++) {
		curr_pin_power = _FSRs_to_pin_powers[i];

		/* If this pin power is unique and nozero (doesn't match the previous
		 * pin's power), then tally it
		 */
		if (curr_pin_power > 0 && curr_pin_power != prev_pin_power) {
			tot_pin_power += curr_pin_power;
			num_nonzero_pins++;
			prev_pin_power = curr_pin_power;
		}
	}

	/* Compute the average pin power */
	avg_pin_power = tot_pin_power / num_nonzero_pins;

	/* Normalize each pin power to the average non-zero pin power */
	for (int i=0; i < _num_FSRs; i++) {
		_FSRs_to_pin_powers[i] /= avg_pin_power;
	}


	log_printf(NORMAL, "Plotting pin powers...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, _FSRs_to_pin_powers);
	plot(bitMap, "pin_powers", _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);

	return;
}

/*
 * Plot the fluxes for each FSR
 */
void Solver::plotFluxes(int iter_num){

	/* create BitMaps for plotting */
	BitMap<int>* bitMapFSR = new BitMap<int>;
	BitMap<float>* bitMap = new BitMap<float>;
	bitMapFSR->pixel_x = _plotter->getBitLengthX();
	bitMapFSR->pixel_y = _plotter->getBitLengthX();
	bitMap->pixel_x = _plotter->getBitLengthX();
	bitMap->pixel_y = _plotter->getBitLengthY();
	initialize(bitMapFSR);
	initialize(bitMap);
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMapFSR->color_type = RANDOM;
	bitMap->color_type = SCALED;

	/* make FSR BitMap */
	_plotter->makeFSRMap(bitMapFSR->pixels);

//	for (int i = 0; i < NUM_ENERGY_GROUPS; i++){
//
//		std::stringstream string;
//		string << "flux" << i + 1 << "group";
//		std::string title_str = string.str();
//
//		log_printf(DEBUG, "Plotting group %d flux...", (i+1));
//		_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, _FSRs_to_fluxes[i]);
//		plot(bitMap, title_str, _plotter->getExtension());
//	}


	std::stringstream string;
	string << "flux_tot_i_" << iter_num;
	std::string title_str = string.str();

	log_printf(DEBUG, "Plotting total flux...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, _FSRs_to_fluxes[NUM_ENERGY_GROUPS]);
	plot(bitMap, title_str, _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);
}


//void Solver::fixedSourceIteration(int max_iterations) {
//
//	Track* track;
//	int num_segments;
//	std::vector<segment*> segments;
//	double* weights;
//	segment* segment;
//	double* polar_fluxes;
//	double* scalar_flux;
//	double* old_scalar_flux;
//	double* sigma_t;
//	FlatSourceRegion* fsr;
//	double fsr_flux[NUM_ENERGY_GROUPS];
//	double* ratios;
//	double delta;
//	double volume;
//	int t, j, k, s, p, e, pe;
//	int num_threads = _num_azim / 2;
//
//#if !STORE_PREFACTORS
//	double sigma_t_l;
//	int index;
//#endif
//
//	log_printf(INFO, "Fixed source iteration with max_iterations = %d and "
//			"# threads = %d", max_iterations, num_threads);
//
//	/* Loop for until converged or max_iterations is reached */
//	for (int i = 0; i < max_iterations; i++) {
//
//		/* Initialize flux in each region to zero */
//		zeroFSRFluxes();
//		zeroMeshCells();
//
//		/* Loop over azimuthal each thread and azimuthal angle*
//		 * If we are using OpenMP then we create a separate thread
//		 * for each pair of reflecting azimuthal angles - angles which
//		 * wrap into cycles on each other */
//		#if USE_OPENMP && STORE_PREFACTORS
//		#pragma omp parallel for num_threads(num_threads) \
//				private(t, k, j, i, s, p, e, pe, track, segments, \
//						num_segments, weights, polar_fluxes,\
//						segment, fsr, ratios, delta, fsr_flux)
//		#elif USE_OPENMP && !STORE_PREFACTORS
//		#pragma omp parallel for num_threads(num_threads) \
//				private(t, k, j, i, s, p, e, pe, track, segments, \
//						num_segments, weights, polar_fluxes,\
//						segment, fsr, ratios, delta, fsr_flux,\
//						sigma_t_l, index)
//		#endif
//		/* Loop over each thread */
//		for (t=0; t < num_threads; t++) {
//
//			/* Loop over the pair of azimuthal angles for this thread */
//			j = t;
//			while (j < _num_azim) {
//
//			/* Loop over all tracks for this azimuthal angles */
//			for (k = 0; k < _num_tracks[j]; k++) {
//
//				/* Initialize local pointers to important data structures */
//				track = &_tracks[j][k];
//				segments = track->getSegments();
//				num_segments = track->getNumSegments();
//				weights = track->getPolarWeights();
//				polar_fluxes = track->getPolarFluxes();
//
//				/* Loop over each segment in forward direction */
//				for (s = 0; s < num_segments; s++) {
//					segment = segments.at(s);
//					fsr = &_flat_source_regions[segment->_region_id];
//					ratios = fsr->getRatios();
//
//					/* Zero out temporary FSR flux array */
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
//						fsr_flux[e] = 0.0;
//
//
//					/* Initialize the polar angle and energy group counter */
//					pe = 0;
//
//#if !STORE_PREFACTORS
//					sigma_t = segment->_material->getSigmaT();
//
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//						fsr_flux[e] = 0;
//						sigma_t_l = sigma_t[e] * segment->_length;
//						sigma_t_l = std::min(sigma_t_l,10.0);
//						index = sigma_t_l / _pre_factor_spacing;
//						index = std::min(index * 2 * NUM_POLAR_ANGLES,
//												_pre_factor_max_index);
//
//						for (p = 0; p < NUM_POLAR_ANGLES; p++){
//							delta = (polar_fluxes[pe] - ratios[e]) *
//							(1 - (_pre_factor_array[index + 2 * p] * sigma_t_l
//							+ _pre_factor_array[index + 2 * p + 1]));
//							fsr_flux[e] += delta * weights[p];
//							polar_fluxes[pe] -= delta;
//							pe++;
//						}
//					}
//
//#else
//					/* Loop over all polar angles and energy groups */
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
//							delta = (polar_fluxes[pe] -ratios[e]) *
//													segment->_prefactors[e][p];
//							fsr_flux[e] += delta * weights[p];
//							polar_fluxes[pe] -= delta;
//							pe++;
//						}
//					}
//
//#endif
//
//#if CMFD_ACCEL
//					/* if segment crosses a surface in fwd direction, tally current/weight */
//					if (segment->_mesh_surface_fwd != NULL){
//
//						/* set polar angle * energy group to 0 */
//						pe = 0;
//
//						/* increment number of crossings and weight */
//						segment->_mesh_surface_fwd->incrementCrossings(1);
//						segment->_mesh_surface_fwd->incrementWeight(track->getSpacing() / fabs(cos(track->getPhi() - segment->_mesh_surface_fwd->getNormal())), j);
//
//						/* loop over energy groups */
//						for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//							/* loop over polar angles */
//							for (p = 0; p < NUM_POLAR_ANGLES; p++){
//
//								/* increment current (polar and azimuthal weighted flux, group, azimuthal angle)*/
//								segment->_mesh_surface_fwd->incrementCurrent(polar_fluxes[pe]*weights[p]/2.0, e, j);
//
//								log_printf(DEBUG, "fwd tallying cell: %i, surface: %i, energy: %i, current: %f", segment->_mesh_surface_fwd->getMeshCell(), segment->_mesh_surface_fwd->getId(), e, polar_fluxes[pe]*weights[p]/2.0);
//
//								pe++;
//							}
//						}
//					}
//#endif
//
//					/* Increment the scalar flux for this FSR */
//					fsr->incrementFlux(fsr_flux);
//				}
//
//				/* Transfer flux to outgoing track */
//				track->getTrackOut()->setPolarFluxes(track->isReflOut(),0, polar_fluxes);
//
//				/* Loop over each segment in reverse direction */
//				for (s = num_segments-1; s > -1; s--) {
//					segment = segments.at(s);
//					fsr = &_flat_source_regions[segment->_region_id];
//					ratios = fsr->getRatios();
//
//					/* Zero out temporary FSR flux array */
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
//						fsr_flux[e] = 0.0;
//
//
//					/* Initialize the polar angle and energy group counter */
//					pe = GRP_TIMES_ANG;
//
//#if !STORE_PREFACTORS
//					sigma_t = segment->_material->getSigmaT();
//
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//						fsr_flux[e] = 0;
//						sigma_t_l = sigma_t[e] * segment->_length;
//						sigma_t_l = std::min(sigma_t_l,10.0);
//						index = sigma_t_l / _pre_factor_spacing;
//						index = std::min(index * 2 * NUM_POLAR_ANGLES,
//												_pre_factor_max_index);
//
//						for (p = 0; p < NUM_POLAR_ANGLES; p++){
//							delta = (polar_fluxes[pe] - ratios[e]) *
//							(1 - (_pre_factor_array[index + 2 * p] * sigma_t_l
//							+ _pre_factor_array[index + 2 * p + 1]));
//							fsr_flux[e] += delta * weights[p];
//							polar_fluxes[pe] -= delta;
//							pe++;
//						}
//					}
//
//#else
//					/* Loop over all polar angles and energy groups */
//					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
//							delta = (polar_fluxes[pe] - ratios[e]) *
//											segment->_prefactors[e][p];
//							fsr_flux[e] += delta * weights[p];
//							polar_fluxes[pe] -= delta;
//							pe++;
//						}
//					}
//#endif
//
//#if CMFD_ACCEL
//					/* if segment crosses a surface in bwd direction, tally current/weight */
//					if (segment->_mesh_surface_bwd != NULL){
//
//						/* set polar angle * energy group to num groups * num angles */
//						pe = GRP_TIMES_ANG;
//
//						/* increment number of crossings and weight */
//						segment->_mesh_surface_bwd->incrementCrossings(1);
//						segment->_mesh_surface_bwd->incrementWeight(track->getSpacing() / fabs(cos(track->getPhi() - segment->_mesh_surface_bwd->getNormal())), j);
//
//						/* loop over energy groups */
//						for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//							/* loop over polar angles */
//							for (p = 0; p < NUM_POLAR_ANGLES; p++){
//
//								/* increment current (polar and azimuthal weighted flux, group, azimuthal angle)*/
//								segment->_mesh_surface_bwd->incrementCurrent(polar_fluxes[pe]*weights[p]/2.0, e, j);
//
//								log_printf(DEBUG, "bwd tallying cell: %i, surface: %i, energy: %i, current: %f", segment->_mesh_surface_bwd->getMeshCell(), segment->_mesh_surface_bwd->getId(), e, polar_fluxes[pe]*weights[p]/2.0);
//
//
//								pe++;
//							}
//						}
//					}
//#endif
//
//					/* Increment the scalar flux for this FSR */
//					fsr->incrementFlux(fsr_flux);
//				}
//
//				/* Transfer flux to incoming track */
//				track->getTrackIn()->setPolarFluxes(track->isReflIn(),GRP_TIMES_ANG, polar_fluxes);
//			}
//
//			/* Update the azimuthal angle index for this thread
//			 * such that the next azimuthal angle is the one that reflects
//			 * out of the current one. If instead this is the 2nd (final)
//			 * angle to be used by this thread, break loop */
//			if (j < num_threads)
//				j = _num_azim - j - 1;
//			else
//				break;
//
//			}
//		}
//
//
//		/* Add in source term and normalize flux to volume for each region */
//		/* Loop over flat source regions, energy groups */
//		#if USE_OPENMP
//		#pragma omp parallel for private(fsr, scalar_flux, ratios, \
//													sigma_t, volume)
//		#endif
//		for (int r = 0; r < _num_FSRs; r++) {
//			fsr = &_flat_source_regions[r];
//			scalar_flux = fsr->getFlux();
//			ratios = fsr->getRatios();
//			sigma_t = fsr->getMaterial()->getSigmaT();
//			volume = fsr->getVolume();
//
//			for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
//				fsr->setFlux(e, scalar_flux[e] / 2.0);
//				fsr->setFlux(e, FOUR_PI * ratios[e] + (scalar_flux[e] /
//												(sigma_t[e] * volume)));
//			}
//		}
//
//
//		/* Check for convergence if max_iterations > 1 */
//		if (max_iterations > 1) {
//			bool converged = true;
//			#if USE_OPENMP
//			#pragma omp parallel for private(fsr, scalar_flux, old_scalar_flux) shared(converged)
//			#endif
//			for (int r = 0; r < _num_FSRs; r++) {
//				fsr = &_flat_source_regions[r];
//				scalar_flux = fsr->getFlux();
//				old_scalar_flux = fsr->getOldFlux();
//
//				for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
//					if (fabs((scalar_flux[e] - old_scalar_flux[e]) /
//							old_scalar_flux[e]) > FLUX_CONVERGENCE_THRESH )
//						converged = false;
//
//					/* Update old scalar flux */
//					old_scalar_flux[e] = scalar_flux[e];
//				}
//			}
//
//			if (converged)
//				return;
//		}
//
//		/* Update the old scalar flux for each region, energy group */
//		#if USE_OPENMP
//		#pragma omp parallel for private(fsr, scalar_flux, old_scalar_flux)
//		#endif
//		for (int r = 0; r < _num_FSRs; r++) {
//			fsr = &_flat_source_regions[r];
//			scalar_flux = fsr->getFlux();
//			old_scalar_flux = fsr->getOldFlux();
//
//			/* Update old scalar flux */
//			for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
//				old_scalar_flux[e] = scalar_flux[e];
//		}
//	}
//
//	if (max_iterations > 1)
//		log_printf(WARNING, "Scalar flux did not converge after %d iterations",
//															max_iterations);
//
//	return;
//}


void Solver::fixedSourceIteration(int max_iterations) {

	Track* track;
	int num_segments;
	std::vector<segment*> segments;
	double* weights;
	segment* segment;
	double* polar_fluxes;
	double* scalar_flux;
	double* old_scalar_flux;
	double* sigma_t;
	FlatSourceRegion* fsr;
	double fsr_flux[NUM_ENERGY_GROUPS];
	double* ratios;
	double delta;
	double volume;
	int t, j, k, s, p, e, pe;
	int num_threads = _num_azim / 2;

#if !STORE_PREFACTORS
	double sigma_t_l;
	int index;
#endif

	log_printf(INFO, "Fixed source iteration with max_iterations = %d and "
			"# threads = %d", max_iterations, num_threads);

	/* Loop for until converged or max_iterations is reached */
	for (int i = 0; i < max_iterations; i++) {

		/* Initialize flux in each region to zero */
		zeroFSRFluxes();
		zeroMeshCells();

		/* Loop over azimuthal each thread and azimuthal angle*
		 * If we are using OpenMP then we create a separate thread
		 * for each pair of reflecting azimuthal angles - angles which
		 * wrap into cycles on each other */
		#if USE_OPENMP && STORE_PREFACTORS
		#pragma omp parallel for num_threads(num_threads) \
				private(t, k, j, i, s, p, e, pe, track, segments, \
						num_segments, weights, polar_fluxes,\
						segment, fsr, ratios, delta, fsr_flux)
		#elif USE_OPENMP && !STORE_PREFACTORS
		#pragma omp parallel for num_threads(num_threads) \
				private(t, k, j, i, s, p, e, pe, track, segments, \
						num_segments, weights, polar_fluxes,\
						segment, fsr, ratios, delta, fsr_flux,\
						sigma_t_l, index)
		#endif
		/* Loop over each thread */
		for (t=0; t < num_threads; t++) {

			/* Loop over the pair of azimuthal angles for this thread */
			j = t;
			while (j < _num_azim) {

			/* Loop over all tracks for this azimuthal angles */
			for (k = 0; k < _num_tracks[j]; k++) {

				/* Initialize local pointers to important data structures */
				track = &_tracks[j][k];
				segments = track->getSegments();
				num_segments = track->getNumSegments();
				weights = track->getPolarWeights();
				polar_fluxes = track->getPolarFluxes();

				/* Loop over each segment in forward direction */
				for (s = 0; s < num_segments; s++) {
					segment = segments.at(s);
					fsr = &_flat_source_regions[segment->_region_id];
					ratios = fsr->getRatios();

					/* Zero out temporary FSR flux array */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
						fsr_flux[e] = 0.0;


					/* Initialize the polar angle and energy group counter */
					pe = 0;

#if !STORE_PREFACTORS
					sigma_t = segment->_material->getSigmaT();

					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

						fsr_flux[e] = 0;
						sigma_t_l = sigma_t[e] * segment->_length;
						sigma_t_l = std::min(sigma_t_l,10.0);
						index = sigma_t_l / _pre_factor_spacing;
						index = std::min(index * 2 * NUM_POLAR_ANGLES,
												_pre_factor_max_index);

						for (p = 0; p < NUM_POLAR_ANGLES; p++){
							delta = (polar_fluxes[pe] - ratios[e]) *
							(1 - (_pre_factor_array[index + 2 * p] * sigma_t_l
							+ _pre_factor_array[index + 2 * p + 1]));
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#else
					/* Loop over all polar angles and energy groups */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
							delta = (polar_fluxes[pe] -ratios[e]) *
													segment->_prefactors[e][p];
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#endif

#if CMFD_ACCEL
					/* if segment crosses a surface in fwd direction, tally current/weight */
					if (segment->_mesh_surface_fwd != NULL){

						/* set polar angle * energy group to 0 */
						pe = 0;

						/* increment number of crossings and weight */
						segment->_mesh_surface_fwd->incrementCrossings(1);
						segment->_mesh_surface_fwd->incrementWeight(track->getSpacing() / fabs(cos(track->getPhi() - segment->_mesh_surface_fwd->getNormal())), j);

						/* loop over energy groups */
						for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

							/* loop over polar angles */
							for (p = 0; p < NUM_POLAR_ANGLES; p++){

								/* increment current (polar and azimuthal weighted flux, group, azimuthal angle)*/
								segment->_mesh_surface_fwd->incrementCurrent(polar_fluxes[pe]*weights[p]/2.0, e, j);

								log_printf(DEBUG, "fwd tallying cell: %i, surface: %i, energy: %i, current: %f", segment->_mesh_surface_fwd->getMeshCell(), segment->_mesh_surface_fwd->getId(), e, polar_fluxes[pe]*weights[p]/2.0);

								pe++;
							}
						}
					}
#endif

					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}

				/* Transfer flux to outgoing track */
				track->getTrackOut()->setNewPolarFluxes(track->isReflOut(),0, polar_fluxes);

				/* Loop over each segment in reverse direction */
				for (s = num_segments-1; s > -1; s--) {
					segment = segments.at(s);
					fsr = &_flat_source_regions[segment->_region_id];
					ratios = fsr->getRatios();

					/* Zero out temporary FSR flux array */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
						fsr_flux[e] = 0.0;


					/* Initialize the polar angle and energy group counter */
					pe = GRP_TIMES_ANG;

#if !STORE_PREFACTORS
					sigma_t = segment->_material->getSigmaT();

					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

						fsr_flux[e] = 0;
						sigma_t_l = sigma_t[e] * segment->_length;
						sigma_t_l = std::min(sigma_t_l,10.0);
						index = sigma_t_l / _pre_factor_spacing;
						index = std::min(index * 2 * NUM_POLAR_ANGLES,
												_pre_factor_max_index);

						for (p = 0; p < NUM_POLAR_ANGLES; p++){
							delta = (polar_fluxes[pe] - ratios[e]) *
							(1 - (_pre_factor_array[index + 2 * p] * sigma_t_l
							+ _pre_factor_array[index + 2 * p + 1]));
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#else
					/* Loop over all polar angles and energy groups */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
							delta = (polar_fluxes[pe] - ratios[e]) *
											segment->_prefactors[e][p];
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}
#endif

#if CMFD_ACCEL
					/* if segment crosses a surface in bwd direction, tally current/weight */
					if (segment->_mesh_surface_bwd != NULL){

						/* set polar angle * energy group to num groups * num angles */
						pe = GRP_TIMES_ANG;

						/* increment number of crossings and weight */
						segment->_mesh_surface_bwd->incrementCrossings(1);
						segment->_mesh_surface_bwd->incrementWeight(track->getSpacing() / fabs(cos(track->getPhi() - segment->_mesh_surface_bwd->getNormal())), j);

						/* loop over energy groups */
						for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

							/* loop over polar angles */
							for (p = 0; p < NUM_POLAR_ANGLES; p++){

								/* increment current (polar and azimuthal weighted flux, group, azimuthal angle)*/
								segment->_mesh_surface_bwd->incrementCurrent(polar_fluxes[pe]*weights[p]/2.0, e, j);

								log_printf(DEBUG, "bwd tallying cell: %i, surface: %i, energy: %i, current: %f", segment->_mesh_surface_bwd->getMeshCell(), segment->_mesh_surface_bwd->getId(), e, polar_fluxes[pe]*weights[p]/2.0);


								pe++;
							}
						}
					}
#endif

					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}

				/* Transfer flux to incoming track */
				track->getTrackIn()->setNewPolarFluxes(track->isReflIn(),GRP_TIMES_ANG, polar_fluxes);
			}

			/* Update the azimuthal angle index for this thread
			 * such that the next azimuthal angle is the one that reflects
			 * out of the current one. If instead this is the 2nd (final)
			 * angle to be used by this thread, break loop */
			if (j < num_threads)
				j = _num_azim - j - 1;
			else
				break;

			}
		}

		/* set track polar fluxes to new polar fluxes */

		/* loop over azimuthal angles */
		for (int azim = 0; azim < _num_azim; azim++){

			/* Loop over all tracks for this azimuthal angles */
			for (int tr = 0; tr < _num_tracks[azim]; tr++) {

				track = &_tracks[azim][tr];
				track->setPolarFluxes(false,0,track->getNewPolarFluxes());
				track->setPolarFluxes(true,GRP_TIMES_ANG,track->getNewPolarFluxes());

			}
		}


		/* Add in source term and normalize flux to volume for each region */
		/* Loop over flat source regions, energy groups */
		#if USE_OPENMP
		#pragma omp parallel for private(fsr, scalar_flux, ratios, \
													sigma_t, volume)
		#endif
		for (int r = 0; r < _num_FSRs; r++) {
			fsr = &_flat_source_regions[r];
			scalar_flux = fsr->getFlux();
			ratios = fsr->getRatios();
			sigma_t = fsr->getMaterial()->getSigmaT();
			volume = fsr->getVolume();

			for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
				fsr->setFlux(e, scalar_flux[e] / 2.0);
				fsr->setFlux(e, FOUR_PI * ratios[e] + (scalar_flux[e] /
												(sigma_t[e] * volume)));
			}
		}


		/* Check for convergence if max_iterations > 1 */
		if (max_iterations > 1) {
			bool converged = true;
			#if USE_OPENMP
			#pragma omp parallel for private(fsr, scalar_flux, old_scalar_flux) shared(converged)
			#endif
			for (int r = 0; r < _num_FSRs; r++) {
				fsr = &_flat_source_regions[r];
				scalar_flux = fsr->getFlux();
				old_scalar_flux = fsr->getOldFlux();

				for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
					if (fabs((scalar_flux[e] - old_scalar_flux[e]) /
							old_scalar_flux[e]) > FLUX_CONVERGENCE_THRESH )
						converged = false;

					/* Update old scalar flux */
					old_scalar_flux[e] = scalar_flux[e];
				}
			}

			if (converged)
				return;
		}

		/* Update the old scalar flux for each region, energy group */
		#if USE_OPENMP
		#pragma omp parallel for private(fsr, scalar_flux, old_scalar_flux)
		#endif
		for (int r = 0; r < _num_FSRs; r++) {
			fsr = &_flat_source_regions[r];
			scalar_flux = fsr->getFlux();
			old_scalar_flux = fsr->getOldFlux();

			/* Update old scalar flux */
			for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
				old_scalar_flux[e] = scalar_flux[e];
		}
	}

	if (max_iterations > 1)
		log_printf(WARNING, "Scalar flux did not converge after %d iterations",
															max_iterations);

	return;
}



/* initialize the source and renormalize the fsr and track fluxes */
void Solver::initializeSource(){

	double scatter_source, fission_source;
	double fis_source_tot, abs_source_tot;
	double renorm_factor, volume;
	double* nu_sigma_f;
	double* sigma_a;
	double* sigma_s;
	double* chi;
	double* scalar_flux;
	double* source;
	FlatSourceRegion* fsr;
	Material* material;
	int start_index, end_index;


	/*********************************************************************
	 * Renormalize scalar and boundary fluxes
	 *********************************************************************/

	/* Initialize fission source to zero */
	fission_source = 0;

	/* Compute total fission source for this region */
	for (int r = 0; r < _num_FSRs; r++) {

		/* Get pointers to important data structures */
		fsr = &_flat_source_regions[r];
		material = fsr->getMaterial();
		nu_sigma_f = material->getNuSigmaF();
		scalar_flux = fsr->getFlux();
		volume = fsr->getVolume();

		start_index = fsr->getMaterial()->getNuSigmaFStart();
		end_index = fsr->getMaterial()->getNuSigmaFEnd();

		for (int e = start_index; e < end_index; e++)
			fission_source += nu_sigma_f[e] * scalar_flux[e] * volume;
	}

	/* Renormalize scalar fluxes in each region */
	renorm_factor = 1.0 / fission_source;

	#if USE_OPENMP
	#pragma omp parallel for
	#endif
	for (int r = 0; r < _num_FSRs; r++)
		_flat_source_regions[r].normalizeFluxes(renorm_factor);

	/* Renormalization angular boundary fluxes for each track */
	#if USE_OPENMP
	#pragma omp parallel for
	#endif
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++)
			_tracks[i][j].normalizeFluxes(renorm_factor);
	}


	/*********************************************************************
	 * Compute the source for each region
	 *********************************************************************/

	abs_source_tot = 0;
	fis_source_tot = 0;

	/* For all regions, find the source */
	for (int r = 0; r < _num_FSRs; r++) {

		fsr = &_flat_source_regions[r];
		material = fsr->getMaterial();

		/* Initialize the fission source to zero for this region */
		fission_source = 0;
		scalar_flux = fsr->getFlux();
		source = fsr->getSource();
		material = fsr->getMaterial();
		nu_sigma_f = material->getNuSigmaF();
		sigma_a    = material->getSigmaA();
		chi = material->getChi();
		sigma_s = material->getSigmaS();

		start_index = material->getNuSigmaFStart();
		end_index = material->getNuSigmaFEnd();

		/* Compute total fission source for current region */
		for (int e = start_index; e < end_index; e++)
			fission_source += scalar_flux[e] * nu_sigma_f[e];

		/* Compute total scattering source for group G */
		for (int G = 0; G < NUM_ENERGY_GROUPS; G++) {
			abs_source_tot  += scalar_flux[G] * sigma_a[G] * fsr->getVolume();
			fis_source_tot  += scalar_flux[G] * nu_sigma_f[G] * fsr->getVolume();
			scatter_source = 0;

			start_index = material->getSigmaSStart(G);
			end_index = material->getSigmaSEnd(G);

			for (int g = start_index; g < end_index; g++)
				scatter_source += sigma_s[G*NUM_ENERGY_GROUPS + g]
				                          * scalar_flux[g];

			/* Set the total source for region r in group G */
			source[G] = ((1.0 / (_old_k_effs.front())) * fission_source *
							chi[G] + scatter_source) * ONE_OVER_FOUR_PI;
		}
	}

	/* Update pre-computed source / sigma_t ratios */
	computeRatios();
}


double Solver::computeKeff(int max_iterations) {

	log_printf(NORMAL, "Computing k_eff...");

	double* source;
	double* old_source;
	FlatSourceRegion* fsr;

	double cmfd_keff = 0.0;

	/* Check that each FSR has at least one segment crossing it */
	checkTrackSpacing();

	/* Initial guess */
	_old_k_effs.push(1.0);

	/* Set scalar flux to unity for each region */
	oneFSRFluxes();
	zeroTrackFluxes();

	/* Set the old source to unity for each Region */
	#if USE_OPENMP
	#pragma omp parallel for private(fsr)
	#endif
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_flat_source_regions[r];

		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			fsr->setOldSource(e, 1.0);
	}

	/* Source iteration loop */
	for (int i = 0; i < max_iterations; i++) {

		log_printf(NORMAL, "Iteration %d: k_eff = %f", i, _k_eff);

		/* initialize the source and renormalize fsr and track fluxes */
		initializeSource();

		/* Iterate on the flux with the new source */
		fixedSourceIteration(1);

#if CMFD_ACCEL
		_geom->getMesh()->splitCorners(_num_azim);

		/* compute the total surface currents */
		_geom->getMesh()->computeTotCurrents(_num_azim);

		/* compute mesh cell XSs */
		computeXS(_geom->getMesh());

		/* compute mesh cell Ds */
		computeDs(_geom->getMesh());

		/* Perform eigenvalue solve to get CMFD keff */
		cmfd_keff = computeCMFDFlux(_geom->getMesh(), _k_eff);

		checkNeutBal(_geom->getMesh(), cmfd_keff);

		if (_plotter->plotFlux() == true){
			/* Load fluxes into FSR to flux map */
			for (int r=0; r < _num_FSRs; r++) {
				double* fluxes = _flat_source_regions[r].getFlux();
				for (int e=0; e < NUM_ENERGY_GROUPS; e++){
					_FSRs_to_fluxes[e][r] = fluxes[e];
					_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] =
						_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] + fluxes[e];
				}
			}
			plotFluxes(i);
			_plotter->plotCMFDflux(_geom->getMesh(), i);
			_plotter->plotXS(_geom->getMesh(), i);
			_plotter->plotDHats(_geom->getMesh(), i);
		}


#endif
		/* Update k_eff */
		if (_update_flux == true){
			updateMOCFlux(_geom->getMesh());
			_k_eff = cmfd_keff;
		}
		else{
			updateKeff();
		}



		/* If k_eff converged, return k_eff */
		if (fabs(_old_k_effs.back() - _k_eff) < KEFF_CONVERG_THRESH){

			/* Plot net current, surface flux, xs, and d_hats for mesh */
			if (_plotter->plotCurrent()){
#if CMFD_ACCEL
				_plotter->plotNetCurrents(_geom->getMesh());
//				_plotter->plotXS(_geom->getMesh());
//				_plotter->plotDHats(_geom->getMesh());
#endif
			}
			else {
				/* Converge the scalar flux spatially within geometry to plot */
				fixedSourceIteration(1000);
			}

			if (_plotter->plotFlux() == true){
				/* Load fluxes into FSR to flux map */
				for (int r=0; r < _num_FSRs; r++) {
					double* fluxes = _flat_source_regions[r].getFlux();
					for (int e=0; e < NUM_ENERGY_GROUPS; e++){
						_FSRs_to_fluxes[e][r] = fluxes[e];
						_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] =
							_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] + fluxes[e];
					}
				}
				plotFluxes(i+1);
			}


			return _k_eff;
		}

		/* If not converged, save old k_eff and pop off old_keff from
		 *  previous iteration */
		_old_k_effs.push(_k_eff);
		if (_old_k_effs.size() == NUM_KEFFS_TRACKED)
			_old_k_effs.pop();

		/* Update sources in each FSR */
		for (int r = 0; r < _num_FSRs; r++) {
			fsr = &_flat_source_regions[r];
			source = fsr->getSource();
			old_source = fsr->getOldSource();

			for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
				old_source[e] = source[e];
		}
	}

	log_printf(WARNING, "Unable to converge the source after %d iterations",
															max_iterations);

	/* Converge the scalar flux spatially within geometry to plot */
	fixedSourceIteration(1000);

	if (_plotter->plotFlux() == true){
		log_printf(NORMAL, "Plotting fluxes...");
		/* Load fluxes into FSR to flux map */
		for (int r=0; r < _num_FSRs; r++) {
			double* fluxes = _flat_source_regions[r].getFlux();
			for (int e=0; e < NUM_ENERGY_GROUPS; e++){
				_FSRs_to_fluxes[e][r] = fluxes[e];
				_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] =
						_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] + fluxes[e];
			}
		}
		plotFluxes(100);
	}

	return _k_eff;
}


/* compute the cross section for all MeshCells in the Mesh */
void Solver::computeXS(Mesh* mesh){

	log_printf(NORMAL, "Computing cmfd cross sections...");

	/* initialize variables */
	double volume, flux, abs, tot, nu_fis, chi;
	double* scat;
	double abs_tally_cell, flux_tally_cell;
	double nu_fis_tally_cell, dif_tally_cell, rxn_tally_cell, vol_tally_cell;
	double scat_tally_cell[NUM_ENERGY_GROUPS];

	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	Material* material;

	/* loop over mesh cells */
	for (int i = 0; i < mesh->getCellWidth() * mesh->getCellHeight(); i++){
		meshCell = mesh->getCells(i);

		/* loop over energy groups */
		std::vector<int>::iterator iter;
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			abs_tally_cell = 0;
			flux_tally_cell = 0;
			nu_fis_tally_cell = 0;
			dif_tally_cell = 0;
			rxn_tally_cell = 0;
			vol_tally_cell = 0;

			for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
				scat_tally_cell[g] = 0;
			}

			/* loop over FSRs in mesh cell */
			for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter){
				fsr = &_flat_source_regions[*iter];

				/* Gets FSR specific data. */
				material = fsr->getMaterial();
				chi = material->getChi()[e];
				volume = fsr->getVolume();
				flux = fsr->getFlux()[e];

				abs = material->getSigmaA()[e];
				tot = material->getSigmaT()[e];
				nu_fis = material->getNuSigmaF()[e];
				scat = material->getSigmaS();

				log_printf(DEBUG, "tallying fsr = %i, group = %i, flux = %f, volume = %f, mesh cell = %i", *iter, e, flux, volume, i);

				/* Tally rxn rates for FSR in cell */
				abs_tally_cell += abs * flux * volume;
				nu_fis_tally_cell += nu_fis * flux * volume;
				flux_tally_cell += flux;
				dif_tally_cell += flux  * volume / (3.0 * tot);
				rxn_tally_cell += flux * volume;
				vol_tally_cell += volume;

				log_printf(DEBUG, "abs: %f, nu_fis: %f, flux: %f, vol: %f, tot: %f", abs, nu_fis, flux, volume, tot);

				for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
					scat_tally_cell[g] += scat[g*NUM_ENERGY_GROUPS + e] * flux * volume;
					log_printf(DEBUG, "e: %i, g: %i, scat: %f", e,g,scat[g*NUM_ENERGY_GROUPS+e]);
				}

				if (chi >= meshCell->getChi()[e]){
					meshCell->setChi(chi,e);
				}

			}

			meshCell->setVolume(vol_tally_cell);
			meshCell->setSigmaA(abs_tally_cell / rxn_tally_cell, e);
			meshCell->setNuSigmaF(nu_fis_tally_cell / rxn_tally_cell, e);
			meshCell->setDiffusivity(dif_tally_cell / rxn_tally_cell, e);
			meshCell->setOldFlux(rxn_tally_cell / vol_tally_cell, e);
			log_printf(DEBUG, "cell: %i, e: %i, vol: %f, sigA: %f, nuSigF: %f, flux: %f, dif: %f", i, e, vol_tally_cell, abs_tally_cell / rxn_tally_cell,
					nu_fis_tally_cell / rxn_tally_cell, rxn_tally_cell / vol_tally_cell, dif_tally_cell / rxn_tally_cell);

			for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
				meshCell->setSigmaS(scat_tally_cell[g] / rxn_tally_cell,e,g);
			}

		}
	}


	/* compute fis and abs rates based on tallied fluxes an XSs */
	double fis_tot = 0;
	double abs_tot = 0;
	for (int i = 0; i < mesh->getCellWidth() * mesh->getCellHeight(); i++){
		meshCell = mesh->getCells(i);
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
			fis_tot += meshCell->getNuSigmaF()[e]*meshCell->getOldFlux()[e]*meshCell->getVolume();
			abs_tot += meshCell->getSigmaA()[e]*meshCell->getOldFlux()[e]*meshCell->getVolume();
		}
	}

	/* print keff based on nu_fis / abs */
	log_printf(NORMAL, "fission rate / abs rate: %f", fis_tot / abs_tot);
}


/* compute the xs for all MeshCells in the Mesh */
void Solver::computeDs(Mesh* mesh){

	/* initialize variables */
	double d, d_next, d_hat, d_tilde, current, flux, flux_next;
	MeshCell* meshCell;
	MeshCell* meshCellNext;

	/* set cell width and height */
	int cell_height = mesh->getCellHeight();
	int cell_width = mesh->getCellWidth();

	/* loop over mesh cells in y direction */
	for (int y = 0; y < cell_height; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cell_width; x++){

			/* get mesh cell */
			meshCell = mesh->getCells(y*cell_width + x);

			for (int e = 0; e < NUM_ENERGY_GROUPS; e++){

				/* get diffusivity and flux for mesh cell */
				d = meshCell->getDiffusivity()[e];
				flux = meshCell->getOldFlux()[e];

				/* compute d_hat and d_tilde for LEFT surface */
				/* if cell on left side, set d_hat and d_tilde to 0 */
				if (x == 0){
					d_hat = 0.0;
					d_tilde = 0.0;
				}
				else{
					/* get mesh cell to left */
					meshCellNext = mesh->getCells(y*cell_width + x - 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* compute d_hat */
					d_hat = 2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next);

					/* get net outward current across surface */
					current = 0.0;
					/* increment current by outward current on next cell's right side */
					current += meshCellNext->getMeshSurfaces(2)->getCurrentTot(e);

					/* decrement current by outward current on left side */
					current -= meshCell->getMeshSurfaces(0)->getCurrentTot(e);

					log_printf(DEBUG, "left side current = %f", current);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current) / (flux_next + flux);
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(0)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(0)->setDTilde(d_tilde, e);

				/* compute d_hat and d_tilde for BOTTOM surface */
				/* if cell on bottom side, set d_hat and d_tilde to 0 */
				if (y == cell_height - 1){
					d_hat = 0.0;
					d_tilde = 0.0;
				}
				else{
					/* get mesh cell below */
					meshCellNext = mesh->getCells((y+1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* compute d_hat */
					d_hat = 2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on bottom side */
					current += meshCell->getMeshSurfaces(1)->getCurrentTot(e);

					/* decrement current by outward current on next cell's top side */
					current -= meshCellNext->getMeshSurfaces(3)->getCurrentTot(e);

					log_printf(DEBUG, "bottom side current = %f", current);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current) / (flux_next + flux);
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(1)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(1)->setDTilde(d_tilde, e);

				/* compute d_hat and d_tilde for RIGHT surface */
				/* if cell on right side, set d_hat and d_tilde to 0 */
				if (x == cell_width - 1){
					d_hat = 0.0;
					d_tilde = 0.0;
				}
				else{

					/* get mesh cell to the right */
					meshCellNext = mesh->getCells(y*cell_width + x + 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* compute d_hat */
					d_hat = 2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on right side */
					current += meshCell->getMeshSurfaces(2)->getCurrentTot(e);

					/* decrement current by outward current on next cell's left side */
					current -= meshCellNext->getMeshSurfaces(0)->getCurrentTot(e);

					log_printf(DEBUG, "right side current = %f", current);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current) / (flux_next + flux);
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(2)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(2)->setDTilde(d_tilde, e);

				/* compute d_hat and d_tilde for TOP surface */
				/* if cell on top side, set d_hat and d_tilde to 0 */
				if (y == 0){
					d_hat = 0.0;
					d_tilde = 0.0;
				}
				else{
					/* get mesh cell above */
					meshCellNext = mesh->getCells((y-1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* compute d_hat */
					d_hat = 2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on next cell's bottom side */
					current += meshCellNext->getMeshSurfaces(1)->getCurrentTot(e);

					/* decrement current by outward current on top side */
					current -= meshCell->getMeshSurfaces(3)->getCurrentTot(e);

					log_printf(DEBUG, "top side current = %f", current);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current) / (flux_next + flux);
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(3)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(3)->setDTilde(d_tilde, e);
			}
		}
	}
}


/*
 * compute the flux in each mesh cell
 */
double Solver::computeCMFDFlux(Mesh* mesh, double keff_moc){

	using namespace arma;

	MeshCell* meshCell;
	double sum_phi_old = 0.0;
	double sum_phi_new = 0.0;
	double keff;
	int ng = NUM_ENERGY_GROUPS;

	int ch = mesh->getCellHeight();
	int cw = mesh->getCellWidth();

	mat A = zeros(ch*cw*ng, ch*cw*ng);
	mat M = zeros(ch*cw*ng, ch*cw*ng);
	vec phi_old = zeros(ch*cw*ng);
	cx_vec eigval;
	cx_mat eigvec;

	/* construct A matrix, M matrix, and flux vector */

	/* loop over mesh cells in y direction */
	for (int y = 0; y < ch; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cw; x++){

			/* loop over energy groups */
			for (int e = 0; e < ng; e++){

				/* get mesh cell */
				meshCell = mesh->getCells(y*cw + x);

				/* diagonal - A */
				A((y*cw + x)*ng + e,(y*cw + x)*ng + e) = meshCell->getSigmaA()[e] * meshCell->getVolume();

				/* scattering out */
				for (int g = 0; g < ng; g++){
					if (e != g){
						log_printf(DEBUG, "scatter out cell: %i, group %i to group %i, sigS: %f", y*cw+x, e, g,meshCell->getSigmaS()[e*ng + g]);
						A((y*cw + x)*ng+e,(y*cw + x)*ng+e) += meshCell->getSigmaS()[e*ng + g] * meshCell->getVolume();
					}
				}

				/* diagonal - M */
				for (int g = 0; g < ng; g++){
					log_printf(DEBUG, "fission from group: %i, to group %i, chi: %f, nuSigF: %f", g, e,meshCell->getChi()[e], meshCell->getNuSigmaF()[g]);
					M((y*cw + x)*ng+e,(y*cw + x)*ng + g) += meshCell->getChi()[e] * meshCell->getNuSigmaF()[g] * meshCell->getVolume();
				}

				/* scattering in */
				for (int g = 0; g < ng; g++){
					if (e != g){
						A((y*cw + x)*ng+e,(y*cw + x)*ng + g) -= meshCell->getSigmaS()[g*ng + e] * meshCell->getVolume();
					}
				}


				/* right */
				if (x != cw - 1){
					/* diagonal */
					A((y*cw + x)*ng + e,(y*cw + x)*ng + e) += (meshCell->getMeshSurfaces(2)->getDHat()[e]      - meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();

					/* transport in */
					A((y*cw + x)*ng + e,(y*cw + x + 1)*ng + e) = - (meshCell->getMeshSurfaces(2)->getDHat()[e] + meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();
				}

				/* left */
				if (x != 0){
					/* diagonal */
					A((y*cw + x)*ng + e,(y*cw + x)*ng + e) += (meshCell->getMeshSurfaces(0)->getDHat()[e]      + meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();

					/* transport in */
					A((y*cw + x)*ng + e,(y*cw + x - 1)*ng + e) = - (meshCell->getMeshSurfaces(0)->getDHat()[e] - meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
				}

				/* below */
				if (y != ch - 1){
					/* diagonal */
					A((y*cw + x)*ng + e,(y*cw + x)*ng + e) += (meshCell->getMeshSurfaces(1)->getDHat()[e]      - meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();

					/* transport in */
					A((y*cw + x)*ng + e,((y+1)*cw + x)*ng + e) = - (meshCell->getMeshSurfaces(1)->getDHat()[e] + meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
				}

				/* above */
				if (y != 0){
					/* diagonal */
					A((y*cw + x)*ng + e,(y*cw + x)*ng + e) += (meshCell->getMeshSurfaces(3)->getDHat()[e]      + meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();

					/* transport in */
					A((y*cw + x)*ng + e,((y-1)*cw + x)*ng + e) = - (meshCell->getMeshSurfaces(3)->getDHat()[e] - meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
				}
			}
		}
	}

	log_printf(NORMAL, "Running CMFD Diffusion Solver using Armadillo eigenvalue solver...");

	/* compute eigenvector and eigenvalues */
	eig_gen(eigval, eigvec, inv(A)*M);

	vec e_val = conv_to<vec>::from(eigval);
	mat e_vec = conv_to<mat>::from(eigvec);
	vec phi_new;

	/* get dominant eigenvalue and corresponding eigenvector */
	for (int i = 0; i < int(e_val.n_rows); i++){
		if (e_val(i) > .01 && e_val(i) < 10){
			keff = double(e_val(i));
			phi_new = e_vec.col(i);
			break;
		}
	}

	log_printf(NORMAL, "CMFD eigenvalue: %f", keff);

	/* normalize flux */
	for (int i = 0; i < cw*ch; i++){
		meshCell = mesh->getCells(i);
		for (int e = 0; e < ng; e++){
			phi_old(i*ng+e) = meshCell->getOldFlux()[e];
		}
	}

//	std::cout << eigval;
//	std::cout << eigvec;
//
//	std::cout << A << std::endl;
//	std::cout << M << std::endl;
//	std::cout << phi_old << std::endl;
//	std::cout << phi_new << std::endl;

	sum_phi_old = sum(phi_old);
	sum_phi_new = sum(phi_new);

	for (int i = 0; i < cw*ch; i++){
		meshCell = mesh->getCells(i);
		for (int e = 0; e < ng; e++){
			meshCell->setNewFlux(phi_new(i*ng + e)*sum_phi_old/sum_phi_new, e);
			phi_new(i*ng + e) = phi_new(i*ng + e)*sum_phi_old/sum_phi_new;
			log_printf(DEBUG, "flux ratio: %f", meshCell->getNewFlux()[e] / meshCell->getOldFlux()[e]);
		}
	}

	log_printf(NORMAL, "CMFD keff: %f", keff);

	return keff;
}

/* update the MOC flux in each FSR */
void Solver::updateMOCFlux(Mesh* mesh){

	/* initialize variables */
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	double old_flux, new_flux, fsr_new_flux;
	double* flux;

	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int ng = NUM_ENERGY_GROUPS;

	/* loop over mesh cells */
	for (int i = 0; i < cw * ch; i++){

		/* get mesh cell and flux values */
		meshCell = mesh->getCells(i);

		for (int e = 0; e < ng; e++){
			old_flux = meshCell->getOldFlux()[e];
			new_flux = meshCell->getNewFlux()[e];
			log_printf(DEBUG, "Updating flux in meshCell: %i, flux ratio: %f", i, new_flux / old_flux);

			/* loop over FRSs in mesh cell */
			std::vector<int>::iterator iter;
			for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter) {
				fsr = &_flat_source_regions[*iter];
				/* get fsr flux */
				flux = fsr->getFlux();
				fsr_new_flux = new_flux / old_flux * flux[e];

				/* set new flux in FSR */
				fsr->setFlux(e, fsr_new_flux);
			}
		}
	}
}


/*
 * check neutron balance in each mesh cell
 */
void Solver::checkNeutBal(Mesh* mesh, double keff){

	log_printf(NORMAL, "Checking neutron balance...");

	/* residual = leakage + absorption - fission */

	double leak, absorb, fis, res, ratio;
	double absorb_tot = 0.0, leak_tot = 0.0, fis_tot = 0.0, ratio_tot = 0.0, res_tot = 0.0, fis_rate = 0.0, fis_rate_tot = 0.0;

	MeshCell* meshCell;
	MeshCell* meshCellNext;

	int cell_height = mesh->getCellHeight();
	int cell_width = mesh->getCellWidth();
	int ng = NUM_ENERGY_GROUPS;

	/* loop over mesh cells in y direction */
	for (int y = 0; y < cell_height; y++){

		// loop over mesh cells in x direction
		for (int x = 0; x < cell_width; x++){

			/* intialize tallies to 0 */
			leak = 0.0;
			absorb = 0.0;
			fis = 0.0;
			res = 0.0;
			ratio = 0.0;

			/* get mesh cell */
			meshCell = mesh->getCells(y*cell_width + x);

			/* cell to right */
			if (x != cell_width - 1){
				meshCellNext = mesh->getCells(y*cell_width + x + 1);
				for (int e = 0; e < ng; e++){
					leak += meshCell->getMeshSurfaces(2)->getCurrentTot(e);
					leak -= meshCellNext->getMeshSurfaces(0)->getCurrentTot(e);
				}
			}

			/* cell to left */
			if (x != 0){
				meshCellNext = mesh->getCells(y*cell_width + x - 1);
				for (int e = 0; e < ng; e++){
					leak += meshCell->getMeshSurfaces(0)->getCurrentTot(e);
					leak -= meshCellNext->getMeshSurfaces(2)->getCurrentTot(e);
				}
			}

			/* above */
			if (y != 0){
				meshCellNext = mesh->getCells((y-1)*cell_width + x);
				for (int e = 0; e < ng; e++){
					leak += meshCell->getMeshSurfaces(3)->getCurrentTot(e);
					leak -= meshCellNext->getMeshSurfaces(1)->getCurrentTot(e);
				}
			}

			/* below */
			if (y != cell_height - 1){
				meshCellNext = mesh->getCells((y+1)*cell_width + x);
				for (int e = 0; e < ng; e++){
					leak += meshCell->getMeshSurfaces(1)->getCurrentTot(e);
					leak -= meshCellNext->getMeshSurfaces(3)->getCurrentTot(e);
				}
			}

			/* compute absorb and fis rates */
			for (int e = 0; e < ng; e ++){
				absorb = meshCell->getSigmaA()[e]*meshCell->getVolume()*meshCell->getNewFlux()[e];
				for (int g = 0; g < ng; g++){
					fis = 1.0 / keff * meshCell->getChi()[e] * meshCell->getNuSigmaF()[g]*meshCell->getVolume()*meshCell->getNewFlux()[g];
					fis_rate = meshCell->getChi()[e] * meshCell->getNuSigmaF()[g]*meshCell->getVolume()*meshCell->getNewFlux()[g];
				}
			}

			/* compute ratio of (fis - abs) / leak and residual */
			ratio = (fis - absorb) / leak;
			res = leak + absorb - fis;

			/* tally total leak, abs, fis, and ratio */
			leak_tot += leak;
			absorb_tot += absorb;
			fis_rate_tot += fis_rate;
			fis_tot += fis;
			ratio_tot += ratio;

			log_printf(NORMAL, "cell: %i, leak: %f, absorb: %f, fis: %f, res: %f, flux: %f, keff: %f, ratio: %f", y*cell_width+x,
					leak, absorb, fis, res, meshCell->getNewFlux()[0], meshCell->getNuSigmaF()[0] / meshCell->getSigmaA()[0], ratio);
		}
	}

	/* compute total residual and average ratio */
	res_tot = leak_tot + absorb_tot - fis_tot;
	ratio_tot = (fis_tot - absorb_tot) / leak_tot;
	log_printf(NORMAL, "Total leak: %f, absorb: %f, fis: %f, res: %f, keff: %f, keff_in: %f, ratio: %f", leak_tot, absorb_tot, fis_tot, res_tot, fis_rate_tot/absorb_tot, keff, ratio_tot);
}




