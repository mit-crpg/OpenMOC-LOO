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
								Plotter* plotter, bool updateFlux) {
	_geom = geom;
	_quad = new Quadrature(TABUCHI);
	_num_FSRs = geom->getNumFSRs();
	_tracks = track_generator->getTracks();
	_num_tracks = track_generator->getNumTracks();
	_num_azim = track_generator->getNumAzim();
	_plotter = plotter;
	_update_flux = updateFlux;
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



void Solver::fixedSourceIteration(int max_iterations, bool cmfd = false) {

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


#if CMFD_ACCEL
//		if (cmfd == true){

			/* zero surface currents */
			Mesh* mesh = _geom->getMesh();
			for (int cell = 0; cell < mesh->getCellHeight()*mesh->getCellWidth(); cell++){
				for (int surface = 0; surface < 8; surface++){
					mesh->getCells(cell)->getMeshSurfaces(surface)->setCrossings(0);
					mesh->getCells(cell)->getMeshSurfaces(surface)->setWeight(0);
					mesh->getCells(cell)->getMeshSurfaces(surface)->setCurrWeight(0);
					for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
						mesh->getCells(cell)->getMeshSurfaces(surface)->setCurrent(0, group);
						mesh->getCells(cell)->getMeshSurfaces(surface)->setFlux(0, group);
					}
				}

			}
//		}
#endif

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
					if (cmfd == true){

						if (segment->_mesh_surface_fwd != NULL){
							pe = 0;
							segment->_mesh_surface_fwd->incrementCrossings(1);
							segment->_mesh_surface_fwd->incrementWeight(track->getAzimuthalWeight());
							segment->_mesh_surface_fwd->incrementCurrWeight(track->getAzimuthalWeight(), track->getPhi());

							for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
								for (p = 0; p < NUM_POLAR_ANGLES; p++){
									/* increment current (polar and azimuthal weighted flux, group)*/
									segment->_mesh_surface_fwd->incrementCurrent(polar_fluxes[pe] * weights[p] * track->getAzimuthalWeight(), track->getPhi(), e);
									segment->_mesh_surface_fwd->incrementFlux(polar_fluxes[pe] * weights[p] * track->getAzimuthalWeight(), e);

									pe++;
								}
							}
						}
					}
#endif


					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}


				/* Transfer flux to outgoing track */
				track->getTrackOut()->setPolarFluxes(track->isReflOut(),
													0, polar_fluxes);

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
					if (cmfd == true){

						if (segment->_mesh_surface_bwd != NULL){
							pe = 0;
							segment->_mesh_surface_bwd->incrementCrossings(1);
							segment->_mesh_surface_bwd->incrementWeight(track->getAzimuthalWeight());
							segment->_mesh_surface_bwd->incrementCurrWeight(track->getAzimuthalWeight(), track->getPhi());

							for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
								for (p = 0; p < NUM_POLAR_ANGLES; p++){
									/* increment current (polar and azimuthal weighted flux, group)*/
									segment->_mesh_surface_bwd->incrementCurrent(polar_fluxes[pe] * weights[p] * track->getAzimuthalWeight(), track->getPhi(), e);
									segment->_mesh_surface_bwd->incrementFlux(polar_fluxes[pe] * weights[p] * track->getAzimuthalWeight(), e);

									pe++;
								}
							}
						}
					}
#endif

					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}

				/* Transfer flux to incoming track */
				track->getTrackIn()->setPolarFluxes(track->isReflIn(),
										GRP_TIMES_ANG, polar_fluxes);
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


double Solver::computeKeff(int max_iterations) {

	double scatter_source, fission_source;
	double renorm_factor, volume;
	double* nu_sigma_f;
	double* sigma_s;
	double* chi;
	double* scalar_flux;
	double* source;
	double* old_source;
	double cmfd_keff;
	FlatSourceRegion* fsr;
	Material* material;
	int start_index, end_index;


	log_printf(NORMAL, "Computing k_eff...");

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

	// Source iteration loop
	for (int i = 0; i < max_iterations; i++) {

		log_printf(NORMAL, "Iteration %d: k_eff = %f", i, _k_eff);

		/*********************************************************************
		 * Renormalize scalar and boundary fluxes
		 *********************************************************************/

		/* Initialize fission source to zero */
		fission_source = 0;

		/* Compute total fission source to zero for this region */
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
			chi = material->getChi();
			sigma_s = material->getSigmaS();

			start_index = material->getNuSigmaFStart();
			end_index = material->getNuSigmaFEnd();

			/* Compute total fission source for current region */
			for (int e = start_index; e < end_index; e++)
				fission_source += scalar_flux[e] * nu_sigma_f[e];

			/* Compute total scattering source for group G */
			for (int G = 0; G < NUM_ENERGY_GROUPS; G++) {
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

		/*********************************************************************
		 * Update flux and check for convergence
		 *********************************************************************/

		/* Update pre-computed source / sigma_t ratios */
		computeRatios();

		/* Iteration the flux with the new source */
		#if CMFD_ACCEL
		fixedSourceIteration(1, true);
		#else
		fixedSourceIteration(1);
		#endif
		/* Update k_eff */
		updateKeff();


		/* If k_eff converged, return k_eff */

		if (fabs(_old_k_effs.back() - _k_eff) < KEFF_CONVERG_THRESH){



			/* Plot net current, surface flux, xs, and d_hats for mesh */
			if (_plotter->plotCurrent()){
				#if CMFD_ACCEL
				fixedSourceIteration(1000, true);
				_geom->getMesh()->splitCorners();
				computeXS(_geom->getMesh());
				computeDs(_geom->getMesh());
				cmfd_keff = computeCMFDFlux(_geom->getMesh());
				log_printf(NORMAL, "final CMFD k_eff: %f", cmfd_keff);
				checkNeutBal(_geom->getMesh(), _k_eff);
				checkNeutBal(_geom->getMesh(), cmfd_keff);
				_plotter->plotNetCurrents(_geom->getMesh());
				_plotter->plotSurfaceFlux(_geom->getMesh());
				_plotter->plotXS(_geom->getMesh());
				_plotter->plotDHats(_geom->getMesh());
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
				plotFluxes();
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

		#if CMFD_ACCEL
		/* Do CMFD solve and update fluxes */
		_geom->getMesh()->splitCorners();
		computeXS(_geom->getMesh());
		computeDs(_geom->getMesh());
		cmfd_keff = computeCMFDFlux(_geom->getMesh());
		checkNeutBal(_geom->getMesh(), _k_eff);
		if (_update_flux == true){
			updateMOCFlux(_geom->getMesh());
		}
		checkNeutBal(_geom->getMesh(), cmfd_keff);
		#endif

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
		plotFluxes();
	}

	return _k_eff;
}


// only plots flux
void Solver::plotFluxes(){

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

	for (int i = 0; i < NUM_ENERGY_GROUPS; i++){

		std::stringstream string;
		string << "flux" << i + 1 << "group";
		std::string title_str = string.str();

		log_printf(NORMAL, "Plotting group %d flux...", (i+1));
		_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, _FSRs_to_fluxes[i]);
		plot(bitMap, title_str, _plotter->getExtension());
	}

	log_printf(NORMAL, "Plotting total flux...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, _FSRs_to_fluxes[NUM_ENERGY_GROUPS]);
	plot(bitMap, "flux_total", _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);
}


/* compute the xs for all MeshCells in the Mesh */
void Solver::computeXS(Mesh* mesh){

	log_printf(NORMAL, "Computing cmfd cross sections...");

	/* tallies for fsr_vol * group_avg_xs and fsr_vol */
	double abs_tally, flux_tally, tot_tally, fis_tally, nu_fis_tally, dif_tally, rxn_tally, vol_tally;
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	Material* material;
	double volume = 0;
	double flux, abs, tot, fis, nu_fis;
	double abs_tally_cell, flux_tally_cell, tot_tally_cell, fis_tally_cell, nu_fis_tally_cell, dif_tally_cell, rxn_tally_cell;


	/* loop over mesh cells */
	for (int i = 0; i < mesh->getCellWidth() * mesh->getCellHeight(); i++){
		meshCell = mesh->getCells(i);
		abs_tally = 0;
		flux_tally = 0;
		tot_tally = 0;
		fis_tally = 0;
		nu_fis_tally = 0;
		dif_tally = 0;
		rxn_tally = 0;

		std::vector<int>::iterator iter;
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			abs_tally_cell = 0;
			flux_tally_cell = 0;
			tot_tally_cell = 0;
			fis_tally_cell = 0;
			nu_fis_tally_cell = 0;
			dif_tally_cell = 0;
			rxn_tally_cell = 0;
			vol_tally = 0;

			for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter){
				fsr = &_flat_source_regions[*iter];
				/* Gets FSR specific data. */
				material = fsr->getMaterial();
				volume = fsr->getVolume();
				flux = fsr->getFlux()[e];
				abs = material->getSigmaA()[e];
				tot = material->getSigmaT()[e];
				fis = material->getSigmaF()[e];
				nu_fis = material->getNuSigmaF()[e];

				/* Tally reaction rates. */
				abs_tally_cell += abs * flux * volume;
				tot_tally_cell += tot * flux * volume;
				fis_tally_cell += fis * flux * volume;
				nu_fis_tally_cell += nu_fis * flux * volume;
				flux_tally_cell += flux;
				dif_tally_cell += flux * volume / (3.0 * tot);
				rxn_tally_cell += flux * volume;
				vol_tally += volume;
			}

			abs_tally += abs_tally_cell;
			tot_tally += tot_tally_cell;
			fis_tally += fis_tally_cell;
			nu_fis_tally += nu_fis_tally_cell;
			dif_tally += dif_tally_cell;
			flux_tally += rxn_tally_cell / vol_tally;
			rxn_tally += rxn_tally_cell;
		}

		meshCell->setSigmaA(abs_tally / flux_tally);
		meshCell->setSigmaT(tot_tally / flux_tally);
		meshCell->setSigmaF(fis_tally / flux_tally);
		meshCell->setSigmaS((meshCell->getSigmaT() - meshCell->getSigmaA()) / flux_tally);
		meshCell->setNuSigmaF(nu_fis_tally / flux_tally);
		meshCell->setDiffusivity(dif_tally / flux_tally);
		meshCell->setOldFlux(flux_tally);
	}
}

/* compute the xs for all MeshCells in the Mesh */
void Solver::computeDs(Mesh* mesh){

	double d, d_next, d_hat, d_tilde, current, flux, flux_next;
	int cell_height, cell_width;
	MeshCell* meshCell;
	MeshCell* meshCellNext;

	cell_height = mesh->getCellHeight();
	cell_width = mesh->getCellWidth();

	/* loop over mesh cells in y direction */
	for (int y = 0; y < cell_height; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cell_width; x++){

			meshCell = mesh->getCells(y*cell_width + x);
			d = meshCell->getDiffusivity();

			/* compute D Hat and D Tilde for left surface */
			if (x == 0){
				d_hat = 0.0;
				d_tilde = 0.0;
			}
			else{
				meshCellNext = mesh->getCells(y*cell_width + x - 1);
				d_next = meshCellNext->getDiffusivity();
				d_hat = 2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next);
				flux = meshCell->getOldFlux();
				flux_next = meshCellNext->getOldFlux();
				current = 0.0;
				for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
					current -= meshCell->getMeshSurfaces(0)->getCurrent(e);
					current += meshCellNext->getMeshSurfaces(2)->getCurrent(e);
				}
				d_tilde = (-d_hat * (flux - flux_next) - current) / (flux_next + flux);
			}
			meshCell->getMeshSurfaces(0)->setDHat(d_hat);
			meshCell->getMeshSurfaces(0)->setDTilde(d_tilde);

			/* compute D Hat for bottom surface */
			if (y == cell_height - 1){
				d_hat = 0.0;
				d_tilde = 0.0;
			}
			else{
				meshCellNext = mesh->getCells((y+1)*cell_width + x);
				d_next = meshCellNext->getDiffusivity();
				d_hat = 2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next);
				flux = meshCell->getOldFlux();
				flux_next = meshCellNext->getOldFlux();
				current = 0.0;
				for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
					current += meshCell->getMeshSurfaces(1)->getCurrent(e);
					current -= meshCellNext->getMeshSurfaces(3)->getCurrent(e);
				}
				d_tilde = (-d_hat * (flux_next - flux) - current) / (flux_next + flux);
			}
			meshCell->getMeshSurfaces(1)->setDHat(d_hat);
			meshCell->getMeshSurfaces(1)->setDTilde(d_tilde);

			/* compute D Hat for right surface */
			if (x == cell_width - 1){
				d_hat = 0.0;
				d_tilde = 0.0;
			}
			else{
				meshCellNext = mesh->getCells(y*cell_width + x + 1);
				d_next = meshCellNext->getDiffusivity();
				d_hat = 2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next);
				flux = meshCell->getOldFlux();
				flux_next = meshCellNext->getOldFlux();
				current = 0.0;
				for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
					current += meshCell->getMeshSurfaces(2)->getCurrent(e);
					current -= meshCellNext->getMeshSurfaces(0)->getCurrent(e);
				}
				d_tilde = (-d_hat * (flux_next - flux) - current) / (flux_next + flux);
			}
			meshCell->getMeshSurfaces(2)->setDHat(d_hat);
			meshCell->getMeshSurfaces(2)->setDTilde(d_tilde);

			/* compute D Hat for top surface */
			if (y == 0){
				d_hat = 0.0;
				d_tilde = 0.0;
			}
			else{
				meshCellNext = mesh->getCells((y-1)*cell_width + x);
				d_next = meshCellNext->getDiffusivity();
				d_hat = 2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next);
				flux = meshCell->getOldFlux();
				flux_next = meshCellNext->getOldFlux();
				current = 0.0;
				for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
					current -= meshCell->getMeshSurfaces(3)->getCurrent(e);
					current += meshCellNext->getMeshSurfaces(1)->getCurrent(e);
				}
				d_tilde = (-d_hat * (flux - flux_next) - current) / (flux_next + flux);
			}
			meshCell->getMeshSurfaces(3)->setDHat(d_hat);
			meshCell->getMeshSurfaces(3)->setDTilde(d_tilde);
		}
	}

}


double Solver::computeCMFDFlux(Mesh* mesh){

	using namespace arma;

	MeshCell* meshCell;
	int cell_height, cell_width;

	cell_height = mesh->getCellHeight();
	cell_width = mesh->getCellWidth();

	mat A = zeros(cell_height*cell_width, cell_height*cell_width);
	mat M = zeros(cell_height*cell_width, cell_height*cell_width);
	colvec phi = ones(cell_height*cell_width, 1);
	colvec sold = zeros(cell_height*cell_width, 1);
	colvec snew = zeros(cell_height*cell_width, 1);
	colvec phi_old = zeros(cell_height*cell_width, 1);

	mat sources = zeros(cell_height*cell_width, 3);

	for (int y = 0; y < cell_height; y++){
		for (int x = 0; x < cell_width; x++){
			phi(y*cell_width+x,0) = mesh->getCells(y*cell_width+x)->getOldFlux();
		}
	}

	/* construct A and M matrix and flux vector */
	// loop over mesh cells in y direction
	for (int y = 0; y < cell_height; y++){

		// loop over mesh cells in x direction
		for (int x = 0; x < cell_width; x++){

			meshCell = mesh->getCells(y*cell_width + x);

			/* diagonal */
			A(y*cell_width + x,y*cell_width + x) = (meshCell->getMeshSurfaces(0)->getDHat() + meshCell->getMeshSurfaces(0)->getDTilde()) * meshCell->getHeight()
					+ (meshCell->getMeshSurfaces(1)->getDHat() - meshCell->getMeshSurfaces(1)->getDTilde()) * meshCell->getWidth()
					+ (meshCell->getMeshSurfaces(2)->getDHat() - meshCell->getMeshSurfaces(2)->getDTilde()) * meshCell->getHeight()
					+ (meshCell->getMeshSurfaces(3)->getDHat() + meshCell->getMeshSurfaces(3)->getDTilde()) * meshCell->getWidth()
					+ meshCell->getSigmaA()*(meshCell->getHeight() * meshCell->getWidth());
			M(y*cell_width + x,y*cell_width + x) = meshCell->getNuSigmaF() * (meshCell->getHeight() * meshCell->getWidth());

			/* upper diagonal */
			if (x != cell_width - 1){
				A(y*cell_width + x,y*cell_width + x + 1) = - (meshCell->getMeshSurfaces(2)->getDHat() + meshCell->getMeshSurfaces(2)->getDTilde()) * meshCell->getHeight();
			}

			/* lower diagonal */
			if (x != 0){
				A(y*cell_width + x,y*cell_width + x - 1) = - (meshCell->getMeshSurfaces(0)->getDHat() - meshCell->getMeshSurfaces(0)->getDTilde()) * meshCell->getHeight();
			}

			/* 2D upper diagonal */
			if (y != cell_height - 1){
				A(y*cell_width + x,(y+1)*cell_width + x) = - (meshCell->getMeshSurfaces(1)->getDHat() + meshCell->getMeshSurfaces(1)->getDTilde()) * meshCell->getWidth();
			}

			/* 2D lower diagonal */
			if (y != 0){
				A(y*cell_width + x,(y-1)*cell_width + x) = - (meshCell->getMeshSurfaces(3)->getDHat() - meshCell->getMeshSurfaces(3)->getDTilde()) * meshCell->getWidth();
			}
		}
	}

	sources.col(0) = A * phi;

	log_printf(NORMAL, "Running CMFD Diffusion Solver using power iteration...");

	int max_outer = 1000;
	double eps = 0.0;
	double keff;
	sold = M * phi;
	double sumold, sumnew;
	sumold = sum(sold);
	sold = sold/(sumold/(cell_width * cell_height));
	phi = inv(A) * sold;
 	for (int y = 0; y < cell_height; y++){
		for (int x = 0; x < cell_width; x++){
			mesh->getCells(y*cell_width+x)->setOldFlux(phi(y*cell_width+x,0));
		}
	}

	sumold = cell_width * cell_height;
	double criteria = 1e-12;
	int iter;

	for (iter = 0; iter < max_outer; iter++){
		phi = inv(A) * sold;
		snew = M * phi;
		sumnew = sum(snew);
		keff = sumnew / sumold;
		eps = 0.0;
		for (int np = 0; np < cell_width * cell_height; np++){
			if (snew(np,0) > 0){
				eps += pow((sold(np,0) * keff) / snew(np,0) - 1, 2);
			}
		}
		eps = pow(eps/(cell_width*cell_height), 0.5);
		sold = snew / (sumnew / (cell_width * cell_height));

		if (iter > 10 && eps < criteria){
			/* update flux in each MeshCell */
			for (int y = 0; y < cell_height; y++){
				for (int x = 0; x < cell_width; x++){
					meshCell = mesh->getCells(y * cell_width + x);
					meshCell->setNewFlux(phi(y*cell_width + x, 0));
					log_printf(NORMAL, "flux ratio: %f", meshCell->getNewFlux() / meshCell->getOldFlux());
				}
			}
			break;
		}
	}

	log_printf(NORMAL, "CMFD keff: %f, iterations: %i", keff, iter);

	sources.col(1) = A * phi;
	sources.col(2) = 1.0 / keff * M * phi;

	std::cout << sources;
	std::cout << A;
	std::cout << M;

	return keff;

}

/* update the MOC flux in each FSR */
void Solver::updateMOCFlux(Mesh* mesh){

	// tallies for fsr_vol * group_avg_xs and fsr_vol
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	double old_flux = 0, new_flux = 0, fsr_new_flux;
	double* flux;

	// loop over mesh cells
	for (int i = 0; i < mesh->getCellWidth() * mesh->getCellHeight(); i++){
		meshCell = mesh->getCells(i);
		old_flux = meshCell->getOldFlux();
		new_flux = meshCell->getNewFlux();
		log_printf(NORMAL, "Updating flux in meshCell: %i, flux ratio: %f", i, new_flux / old_flux);

		std::vector<int>::iterator iter;
		for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter) {
			fsr = &_flat_source_regions[*iter];
			flux = fsr->getFlux();

			for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
				fsr_new_flux = new_flux / old_flux * flux[e];
				fsr->setFlux(e, fsr_new_flux);

			}
		}
	}

}


/* check to make sure the neutron balance equation
 * is satisfied with the currents from the CMFD solver
 */
void Solver::checkNeutBal(Mesh* mesh, double keff){

	log_printf(NORMAL, "Checking neutron balance...");

	/* residual = leakage + absorption - fission */

	double leak, absorb, fis, res;

	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int cell_height, cell_width;

	cell_height = mesh->getCellHeight();
	cell_width = mesh->getCellWidth();

	/* construct A and M matrix and flux vector */
	// loop over mesh cells in y direction
	for (int y = 0; y < cell_height; y++){

		// loop over mesh cells in x direction
		for (int x = 0; x < cell_width; x++){
			meshCell = mesh->getCells(y*cell_width + x);
			leak = 0.0;
			absorb = 0.0;
			fis = 0.0;
			res = 0.0;

			for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
				leak += meshCell->getMeshSurfaces(0)->getCurrent(group)*meshCell->getHeight();
				leak += meshCell->getMeshSurfaces(1)->getCurrent(group)*meshCell->getWidth();
				leak += meshCell->getMeshSurfaces(2)->getCurrent(group)*meshCell->getHeight();
				leak += meshCell->getMeshSurfaces(3)->getCurrent(group)*meshCell->getWidth();
			}

			/* left */
			if (x != 0){
				meshCellNext = mesh->getCells(y*cell_width + x - 1);
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCellNext->getMeshSurfaces(2)->getCurrent(group)*meshCell->getHeight();
				}
			}
			else{
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCell->getMeshSurfaces(0)->getCurrent(group)*meshCell->getHeight();
				}
			}

			/* right */
			if (x != cell_width - 1){
				meshCellNext = mesh->getCells(y*cell_width + x + 1);
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCellNext->getMeshSurfaces(0)->getCurrent(group)*meshCell->getHeight();
				}
			}
			else{
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCell->getMeshSurfaces(2)->getCurrent(group)*meshCell->getHeight();
				}
			}

			/* above */
			if (y != 0){
				meshCellNext = mesh->getCells((y-1)*cell_width + x);
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCellNext->getMeshSurfaces(1)->getCurrent(group)*meshCell->getWidth();
				}
			}
			else{
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCell->getMeshSurfaces(3)->getCurrent(group)*meshCell->getWidth();
				}
			}

			/* below */
			if (y != cell_height - 1){
				meshCellNext = mesh->getCells((y+1)*cell_width + x);
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCellNext->getMeshSurfaces(3)->getCurrent(group)*meshCell->getWidth();
				}
			}
			else{
				for (int group = 0; group < NUM_ENERGY_GROUPS; group++){
					leak -= meshCell->getMeshSurfaces(1)->getCurrent(group)*meshCell->getWidth();
				}
			}

			absorb = meshCell->getSigmaA()*meshCell->getWidth()*meshCell->getHeight()*meshCell->getOldFlux();
			fis = 1.0 / keff * meshCell->getNuSigmaF()*meshCell->getWidth()*meshCell->getHeight()*meshCell->getOldFlux();

			res = leak + absorb - fis;

			log_printf(NORMAL, "cell: %i, leak: %f, absorb: %f, fis: %f, res: %f, flux: %f, keff: %f", y*cell_width+x,
					leak, absorb, fis, res, meshCell->getOldFlux(), meshCell->getNuSigmaF() / meshCell->getSigmaA());

		}
	}





}










