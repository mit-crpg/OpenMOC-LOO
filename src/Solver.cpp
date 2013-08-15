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
               Plotter* plotter, Cmfd* cmfd, Options* opts)
{
    _geom = geom;
    _quad = new Quadrature(TABUCHI);
    _num_FSRs = geom->getNumFSRs();
    _tracks = track_generator->getTracks();
    _num_tracks = track_generator->getNumTracks();
    _num_azim = track_generator->getNumAzim();
    _plotter = plotter;
    _cmfd = cmfd;
    _num_crn = 0;

    /* Options */
    _update_keff = opts->updateKeff();
    _l2_norm_conv_thresh = opts->getL2NormConvThresh();
    _moc_conv_thresh = opts->getMOCConvThresh();
    _compute_powers = opts->computePinPowers();
    _run_cmfd = opts->getCmfd();
    _run_loo = opts->getLoo();
    _run_loo1 = opts->getLoo1();
    _run_loo2 = opts->getLoo2();
    _diffusion = opts->getDiffusion();
    _acc_after_MOC_converge = opts->getAccAfterMOCConverge();
    _k_eff = opts->getKGuess();
    _geometry_file = opts->getGeometryFile();
    _damp_factor = opts->getDampFactor();
    _track_spacing = opts->getTrackSpacing();
    _boundary_iteration = opts->getBoundaryIteration();
    _update_boundary = opts->getUpdateBoundary();
    
    _cmfd_k = _k_eff;
    _loo_k = _k_eff;


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

    /* Initializes FSRs: store id, cell id, material, volume */
    initializeFSRs();

    _total_vol = 0;
    for (int r = 0; r < _num_FSRs; r++) 
        _total_vol += _flat_source_regions[r].getVolume();

    /* Gives cmfd a pointer to the FSRs */
    if (_run_cmfd || _run_loo || _acc_after_MOC_converge) 
    {
        _cw = _geom->getMesh()->getCellWidth();
        _ch = _geom->getMesh()->getCellHeight();
        _cmfd->setFSRs(_flat_source_regions);
        _cmfd->setTracks(_tracks);
        initializeWeights();

        if (_update_boundary)
        {
            log_printf(NORMAL, "Acceleration is on with %d boundary iteration,"
                       " %f damping, and update boundary flux" , 
                       _boundary_iteration, _damp_factor);
        }
        else
        {
            log_printf(NORMAL, "Acceleration is on with %d boundary iteration,"
                       " %f damping, and no update boundary flux" , 
                       _boundary_iteration, _damp_factor);
        }  
    }
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
    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++) 
        {
            curr_track = &_tracks[i][j];
            azim_weight = curr_track->getAzimuthalWeight();

            for (int p = 0; p < NUM_POLAR_ANGLES; p++)
            {
                curr_track->setPolarWeight
                    (p, azim_weight * _quad->getMultiple(p) * FOUR_PI);
            }
        }
    }

    log_printf(DEBUG, "0 azi angle 0 track 0 polar angle omega_m = %.10f", 
               (&_tracks[0][0])->getPolarWeights()[0]);


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
                        curr_seg->_prefactors[e][p] = 
                            computePreFactor(curr_seg, e, p);
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
    int num_array_values = 10 * sqrt(1 / (8 * _moc_conv_thresh));
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
void Solver::initializeFSRs() 
{
    log_printf(NORMAL, "Initializing FSRs...");

    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geom->getUniverse(0);
    Track* track;
    segment* seg;
    FlatSourceRegion* fsr;


    /* Set each FSR's volume by accumulating the total length of all
       tracks inside the FSR. Loop over azimuthal angle, track and segment */
    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++) 
        {
            track = &_tracks[i][j];

            for (int s = 0; s < track->getNumSegments(); s++) 
            {
                seg = track->getSegment(s);
                fsr =&_flat_source_regions[seg->_region_id];
                fsr->incrementVolume(seg->_length * 
                                     track->getAzimuthalWeight());
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
 * Initialize each track's incoming and outgoing polar fluxes
 */
void Solver::initializeTrackFluxes(double flux) {

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
                polar_fluxes[i] = flux;
        }
    }
}


/**
 * Set the scalar flux for each energy group inside each FSR
 * to unity
 */
void Solver::oneFSRFluxOldSource() 
{
    log_printf(INFO, "Setting all FSR scalar fluxes to unity...");
    FlatSourceRegion* fsr;

    /* Loop over all FSRs and energy groups */
#if USE_OPENMP
#pragma omp parallel for private(fsr)
#endif
    for (int r = 0; r < _num_FSRs; r++) 
    {
        fsr = &_flat_source_regions[r];
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            fsr->setFlux(e, 1.0);
            fsr->setOldSource(e, 1.0);
        }
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

    log_printf(INFO, "Setting all mesh cell fluxes and currents to zero...");

    /* get mesh */
    Mesh* mesh = _geom->getMesh();
    MeshCell* meshCell;

    /* loop over mesh cells */
    for (int i = 0; i < mesh->getCellHeight() * mesh->getCellWidth(); i++)
    {
        meshCell = mesh->getCells(i);

        /* loop over mesh surfaces in mesh cell */
        for (int surface = 0; surface < 8; surface++)
        {
            /* loop over energy groups */
            for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
            {
                /* set current to zero */
                meshCell->getMeshSurfaces(surface)->setCurrent(0, e);

                /* set quad currents to zero */
                for (int j = 0; j < 2; j++)
                {
                    meshCell->getMeshSurfaces(surface)->setQuadCurrent(0, e, j);
                }
                
            }
        }
    }
    return;
}


void Solver::zeroLeakage(){

    for (int s = 1; s < 5; s++)
    {
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
            _geom->getSurface(s)->setLeakage(0.0,e);
    }
}

/**
 * Compute k_eff from the new and old source and the value of k_eff from
 * the previous iteration
 */
double Solver::computeKeff(int iteration) 
{
    double tot_abs = 0.0, tot_fission = 0.0, leakage = 0.0;
    double abs = 0, fission = 0, vol = 0;
    double k; 
    double *sigma_a, *nu_sigma_f, *flux;
    Material* material;
    FlatSourceRegion* fsr;

#if USE_OPENMP
#pragma omp parallel shared(tot_abs, tot_fission)
    {
#pragma omp for private(fsr, material, sigma_a, nu_sigma_f, flux, abs, 
	fission)
#endif
	for (int r = 0; r < _num_FSRs; r++) 
	{
            abs = 0;
            fission = 0;
            fsr = &_flat_source_regions[r];
            material = fsr->getMaterial();
            sigma_a = material->getSigmaA();

            nu_sigma_f = material->getNuSigmaF();
            flux = fsr->getFlux();
            vol = fsr->getVolume();

            for (int e = 0; e < NUM_ENERGY_GROUPS; e++) 
            {
                  if (sigma_a[e] * flux[e] * vol < 0.0)
                  {
                      log_printf(DEBUG, "  FSR %d energy %d has abs xs %f, "
                                 " flux %f, vol %f", 
                                 r, e, sigma_a[e], flux[e], vol);
                  }
	        abs += sigma_a[e] * flux[e] * vol;
                fission += nu_sigma_f[e] * flux[e] * vol;
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

for (int s = 1; s < 5; s++)
{
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        leakage += _geom->getSurface(s)->getLeakage()[e];
}

if (leakage < 0.0)
    log_printf(WARNING, 
               "MOC leakage  = %f should be non-negative", leakage);
if (tot_fission < 0.0)
    log_printf(WARNING, "MOC total fission = %f should be positive", 
               tot_fission);
if (tot_abs < 0.0)
    log_printf(WARNING, "MOC total abs = %f should be positive", tot_abs);

k = tot_fission / (tot_abs + leakage);
//_geom->getMesh()->setKeffMOC(_k_eff, iteration);
return k;
}


/**
 * Update FSR's scalar fluxes, normalize them and update the source
 */
void Solver::updateFlux(int moc_iter) 
{
    _cmfd->updateMOCFlux(moc_iter);

    if (_update_boundary)
    {
        if (_run_loo && (!(_diffusion && (moc_iter == 0))))
        {
            log_printf(DEBUG, " iter %d prolongation: update by quadrature",
                       moc_iter);
            if (_run_loo1)
                updateBoundaryFluxByQuadrature();
            else
            {
                updateBoundaryFluxByQuadrature();
                //_cmfd->updateBoundaryFluxBySrc(moc_iter);
                //_cmfd->updateBoundaryFluxByNetCurrent(moc_iter);
                //_cmfd->updateBoundaryFluxByHalfSpace(moc_iter);
            }
        }
        else if ((_diffusion) && (moc_iter == 0))
        {
            _cmfd->updateBoundaryFlux(moc_iter);
        }
        else
        {
            log_printf(DEBUG, " iter %d prolongation: update by partial",
                       moc_iter);
            _cmfd->updateBoundaryFluxByHalfSpace(moc_iter);
        }
    }
    normalizeFlux();
    return;
}

void Solver::updateBoundaryFluxByQuadrature()
{
    Track *track;
    segment *seg;
    MeshSurface *meshSurface, **meshSurfaces;
    double phi, factor;
    int ind, num_segments, pe, surf_id, num_updated = 0;

    meshSurfaces = _geom->getMesh()->getSurfaces();

    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++)
        {
            track = &_tracks[i][j];
            num_segments = track->getNumSegments();
            phi = track->getPhi();
            if (phi < PI / 2.0)
                ind = 1;
            else
                ind = 0;

            /* Forward direction is 0, backward is 1 */
            for (int dir = 0; dir < 2; dir++)
            {
                seg = track->getSegment(dir * (num_segments - 1));

                if (dir == 0)
                    surf_id = seg->_mesh_surface_bwd; 
                else
                    surf_id = seg->_mesh_surface_fwd;

                if ((surf_id % 8) < 4){}
                else if ((surf_id % 8) < 6)
                    surf_id = (surf_id / 8) * 8 + 1;
                else if ((surf_id % 8) < 7)
                    surf_id -= 4;
                else
                    surf_id -= 7;

                meshSurface = meshSurfaces[surf_id];
                pe = dir * GRP_TIMES_ANG;
                for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
                {
                    if (meshSurface->getOldQuadFlux(e, ind) > -100)
                    {
                        factor = meshSurface->getQuadFlux(e, ind) 
                            // / meshSurface->getCurrent(e);
                            / meshSurface->getOldQuadFlux(e, ind);
                        log_printf(DEBUG, "factor = %.10f", factor);
                        for (int p = 0; p < NUM_POLAR_ANGLES; p++)
                        {
                            track->updatePolarFluxes(pe, factor);
                            pe++;
                            num_updated++;
                        }	
                    }	
                    else
                    {
                        if (e == 0)
                        {
                            log_printf(WARNING, 
                                       "e %d i %d j %d (total %d) %f -> %f"
                                       " forward phi %f",
                                       e, i, j, _num_tracks[i],
                                       meshSurface->getOldQuadFlux(e, ind),
                                       meshSurface->getQuadFlux(e, ind), phi);
                        }		
                    }
                }
            }
        }
    }

    log_printf(DEBUG, "updated boundary flux by quadrature per energy: %f", 
               num_updated / (double) NUM_ENERGY_GROUPS);
    return;
}

void Solver::printToScreen(int moc_iter, double eps)
{
    /* Prints & Update keff for MOC sweep */
    if ((moc_iter == 0) && _diffusion)
    {
        printf("Iter %d, MOC k^(m+1) = %.10f, Diffusion k = %.10f,"
               " MOC k^(m+1/2) = %.10f" 
               " FS eps = %.4e,  k eps = %.4e, #Diffusion = %d\n", 
               moc_iter, _k_eff, _cmfd_k, _k_half, eps,
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());

        if (_update_keff)
            _k_eff = _cmfd_k;
    }
    else if ((_run_cmfd) && !(_acc_after_MOC_converge))
    {
        printf("Iter %d, MOC k^(m+1) = %.10f, CMFD k = %.10f,"
               " MOC k^(m+1/2) = %.10f" 
               " FS eps = %.4e,  k eps = %.4e, #CMFD = %d\n", 
               moc_iter, _k_eff, _cmfd_k, _k_half, eps,
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());

        if (_update_keff)
            _k_eff = _cmfd_k;
    }
    else if ((_run_loo) && !(_acc_after_MOC_converge))
    {
#if 1
        printf("Iter %d, MOC k^(m+1) = %.10f, LOO k = %.10f,"
               " MOC k^(m+1/2) = %.10f"
               " FS eps = %.4e, k eps = %.4e, #LOO = %d\n", 
               moc_iter, _k_eff, _loo_k, _k_half, eps,
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());
#else
        log_printf(NORMAL, "Iter %d, MOC k^(m+1) = %.10f, LOO k = %.10f,"
                   " MOC k^(m+1/2) = %.10f"
                   " FS eps = %.4e, k eps = %.4e, #LOO = %d", 
                   moc_iter, _k_eff, _loo_k, _k_half,
                   eps, (_old_k_effs.back() - _old_k_effs.front()) / 
                   _old_k_effs.back(),
                   _cmfd->getNumIterToConv());
#endif
        if (_update_keff)
            _k_eff = _loo_k;
    }
    else
    {
        printf("Iter %d, MOC k = %.10f, FS eps = %.4e, k eps = %.4e\n", 
               moc_iter, _k_eff, eps, 
               (_old_k_effs.back() - _old_k_effs.front()) 
               / _old_k_effs.back());
    }

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
            log_printf(ERROR, 
                       "No tracks were tallied inside FSR id = %d inside "
                       " cell id = %d. This is not going to run because"
                       " the as-tracked vol of this FSR is zero."
                       " Please reduce your track spacing,"
                       " increase the number of azimuthal angles, or increase"
                       " the size of the flat source regions", 
                       i, cell->getId());
        }
    }

    delete [] FSR_segment_tallies;
}

/* Print out the four boundary conditions around the whole geometry in the 
 * order of left, bottom, right, top. 
 */
void Solver::checkBoundary()
{
    std::string bc[4];
    for (int s = 0; s < 4; s++)
    {
        if (_geom->getMesh()->getBoundary(s) == REFLECTIVE)
            bc[s] = "reflective";
        else
            bc[s] = "vacuum";
    }
    log_printf(NORMAL, "Boundary conditions: %s %s %s %s", 
               bc[0].c_str(), bc[1].c_str(), bc[2].c_str(), bc[3].c_str()); 

    return;
}

/**
 * Plot the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates
 */
void Solver::plotPinPowers() {

    log_printf(INFO, "Computing pin powers...");

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
    _plotter->copyFSRMap(bitMapFSR->pixels);

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

void Solver::storeFsrFluxPower() {

    log_printf(INFO, "Computing pin powers...");
    FlatSourceRegion* fsr;
    double *fluxes;

    /* Loop over all FSRs and compute the fission rate*/
    for (int i = 0; i < _num_FSRs; i++) 
    {
        fsr = &_flat_source_regions[i];
        _FSRs_to_powers[i] = fsr->computeFissionRate();
        fluxes = _flat_source_regions[i].getFlux();
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
            _FSRs_to_fluxes[e][i] = fluxes[e];
    }

    return;
}

/*
 * Plot the fluxes for each FSR
 */
void Solver::plotFluxes(int moc_iter)
{
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
    _plotter->copyFSRMap(bitMapFSR->pixels);

    std::stringstream string;
    std::string title_str;
    for (int i = 0; i < NUM_ENERGY_GROUPS; i++)
    {
        string.str("");
        string << "flux" << i + 1 << "group";
        title_str = string.str();
        _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                _FSRs_to_fluxes[i]);
        plot(bitMap, title_str, _plotter->getExtension());
    }

    /* delete bitMaps */
    deleteBitMap(bitMapFSR);
    deleteBitMap(bitMap);
}

void Solver::tallyLooCurrent(Track *track, segment *segment, 
                             MeshSurface **meshSurfaces, int direction)
{
    int index = 0, pe = 0, pe_initial = 0, p, e, surfID;
    MeshSurface *meshSurface;
    double *weights, *polar_fluxes;

    /* Get the ID of the surface that the segment ends on (forward), starts
     * on (backwards)*/
    if (direction == 1)
        surfID = segment->_mesh_surface_fwd;
    else 
    {
        surfID = segment->_mesh_surface_bwd;
        pe_initial = GRP_TIMES_ANG;
    }

    if (surfID != -1)
    {
        /* notice this is polar flux weights, more than just polar weights */
        weights = track->getPolarWeights();
        polar_fluxes = track->getPolarFluxes();

        /* Defines index */
        if (track->getPhi() > PI / 2.0)
            index = 1;
        
        /* Obtains the surface that the segment crosses */
        meshSurface = meshSurfaces[surfID];
        log_printf(DEBUG, " phi = %f, index = %d, surface ID = %d",
                   track->getPhi(), index, 
                   meshSurface->getId());

        pe = pe_initial;
        
        for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
        {
            for (p = 0; p < NUM_POLAR_ANGLES; p++)
            {
                meshSurface->incrementQuadCurrent(polar_fluxes[pe] * weights[p] 
                                                  / 2.0, e, index);
                pe++;
            }
        }
    }
    return;
}

void Solver::tallyCmfdCurrent(Track *track, segment *segment, 
                              MeshSurface **meshSurfaces, int direction)
{
    MeshSurface *meshSurface;
    int pe = 0, pe_initial = 0, p, e, surfID;
    /* polar weights should be azi weight x sin(theta) x polar weight x 4 pi, 
     * where azi weight takes into accout the spacing between tracks */
    double *weights, *polar_fluxes;
    double currents[NUM_ENERGY_GROUPS];
	
    if (direction == 1) 
        surfID = segment->_mesh_surface_fwd;
    else
    {
        surfID = segment->_mesh_surface_bwd;
        pe_initial = GRP_TIMES_ANG;
    }

    if (surfID != -1)
    {

        weights = track->getPolarWeights();
        polar_fluxes = track->getPolarFluxes();

        meshSurface = meshSurfaces[surfID];
        pe = pe_initial;
        for (e = 0; e < NUM_ENERGY_GROUPS; e++)
            currents[e] = 0.0;

        /* loop over energy groups */
        for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
        {
            /* loop over polar angles */
            for (p = 0; p < NUM_POLAR_ANGLES; p++)
            {
                /* increment current (polar flux times polar weights); */
                /* The cos theta needed for getting th direction to normal
                 * is cancelled out with the 1/cos theta needed to get the 
                 * track spacing now on the surface. The 1/2.0 takes into 
                 * account half space. */
                currents[e] += polar_fluxes[pe] * weights[p]/2.0;
                //meshSurface->incrementTotalWt(weights[p] / 2.0, 0);
                pe++;
            }
        }
        meshSurface->incrementCurrent(currents);
    }
    return;
}

void Solver::tallyLooWeight(Track *track, segment *segment, 
                             MeshSurface **meshSurfaces, int direction)
{
    int index = 0, opposite_index = 1, p, surfID;
    MeshSurface *meshSurface;
    double cosTheta, sinTheta, wt, wt2;

    /* Get the ID of the surface that the segment ends on (forward), starts
     * on (backwards)*/
    if (direction == 1)
        surfID = segment->_mesh_surface_fwd;
    else 
        surfID = segment->_mesh_surface_bwd;

    if (surfID != -1)
    {
        /* notice this is polar flux weights, more than just polar weights */
        double *weights = track->getPolarWeights();
        double *sinThetaP = _quad->getSinThetas();

        /* Obtains the surface that the segment crosses */
        meshSurface = meshSurfaces[surfID];

        /* Defines index */
        if (track->getPhi() > PI / 2.0)
        {
            index = 1;
            opposite_index = 0;
        }

        /* Cell ID */
        int i = surfID / 8;
        int y = i / _cw;
        int x = i % _cw;
        int s = surfID % 8;

        cosTheta = fabs(cos(track->getPhi()));
        sinTheta = fabs(sin(track->getPhi()));

        for (p = 0; p < NUM_POLAR_ANGLES; p++)
        {
            wt = 0.5 * weights[p] / TWO_PI / sinThetaP[p];
            wt2 = wt;

            if ((s  == 0) || (s == 2))
            {
                meshSurface->incrementTotalWt(wt2 / cosTheta, index);
                meshSurface->incrementTotalWt(wt, index + 3);
            }
            else if (s < 4)
            {
                meshSurface->incrementTotalWt(wt2 / sinTheta, index);
                meshSurface->incrementTotalWt(wt, index + 3);
            }
            else
            {
                _num_crn += 1;
                wt *= 0.5;
                wt2 *= 0.5;

                if (s < 5)
                {
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt2 / cosTheta, index);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt2 / sinTheta, index);
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt, index + 3);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt, index + 3);

                    if ((x > 0))// && (y < _ch - 1))
                    {
                        meshSurfaces[(i - 1) * 8 + 1]
                            ->incrementTotalWt(wt2 / sinTheta, index);
                        meshSurfaces[(i - 1) * 8 + 1]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 1]
                            ->incrementTotalWt(wt2 / sinTheta, opposite_index);
                        meshSurfaces[i * 8 + 1]
                            ->incrementTotalWt(wt, opposite_index + 3);

                    }

                    if (y < _ch -1) 
                    {
                        meshSurfaces[(i + _cw) * 8 + 0]
                            ->incrementTotalWt(wt2 / cosTheta, index);
                        meshSurfaces[(i + _cw) * 8 + 0]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 0]
                            ->incrementTotalWt(wt2 / cosTheta, opposite_index);
                        meshSurfaces[i * 8 + 0]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }

                }
                else if (s < 6)
                {
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt2 / sinTheta, index);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt2 / cosTheta, index);

                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt, index + 3);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt, index + 3);

                    if (x < _cw - 1)
                    {
                        meshSurfaces[(i + 1) * 8 + 1]
                            ->incrementTotalWt(wt2 / sinTheta, index);
                        meshSurfaces[(i + 1) * 8 + 1]
                            ->incrementTotalWt(wt, index + 3);
                    } 
                    else
                    {
                        meshSurfaces[i * 8 + 1]
                            ->incrementTotalWt(wt2 / sinTheta, opposite_index);
                        meshSurfaces[i * 8 + 1]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    } 
                    
                    if (y < _ch -1) 
                    {
                        meshSurfaces[(i + _cw) * 8 + 2]
                            ->incrementTotalWt(wt2 / cosTheta, index);
                        meshSurfaces[(i + _cw) * 8 + 2]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 2]
                            ->incrementTotalWt(wt2 / cosTheta, opposite_index);
                        meshSurfaces[i * 8 + 2]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }                      
                }
                else if (s < 7)
                {
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt2 / cosTheta, index);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt2 / sinTheta, index);

                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt, index + 3);
                    meshSurfaces[surfID - 3]
                        ->incrementTotalWt(wt, index + 3);

                    if (x < _cw - 1)
                    {
                        meshSurfaces[(i + 1) * 8 + 3]
                            ->incrementTotalWt(wt2 / sinTheta, index);
                        meshSurfaces[(i + 1) * 8 + 3]
                            ->incrementTotalWt(wt, index + 3);
                    } 
                    else
                    {
                        meshSurfaces[i * 8 + 3]
                            ->incrementTotalWt(wt2 / sinTheta, opposite_index);
                        meshSurfaces[i * 8 + 3]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }                      

                    if (y > 0) 
                    {
                        meshSurfaces[(i - _cw) * 8 + 2]
                            ->incrementTotalWt(wt2 / cosTheta, index);
                        meshSurfaces[(i - _cw) * 8 + 2]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 2]
                            ->incrementTotalWt(wt2 / cosTheta, opposite_index);
                        meshSurfaces[i * 8 + 2]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }
                }
                else
                {
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt2 / sinTheta, index);
                    meshSurfaces[surfID - 7]
                        ->incrementTotalWt(wt2 / cosTheta, index);
                    meshSurfaces[surfID - 4]
                        ->incrementTotalWt(wt, index + 3);
                    meshSurfaces[surfID - 7]
                        ->incrementTotalWt(wt, index + 3);


                    if ((x > 0))// && (y > 0))
                    {
                        meshSurfaces[(i - 1) * 8 + 3]
                            ->incrementTotalWt(wt2 / sinTheta, index);
                        meshSurfaces[(i - 1) * 8 + 3]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 3]
                            ->incrementTotalWt(wt2 / sinTheta, opposite_index);
                        meshSurfaces[i * 8 + 3]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }

                    if (y > 0) 
                    {
                        meshSurfaces[(i - _cw) * 8 + 0]
                            ->incrementTotalWt(wt2 / cosTheta, index);
                        meshSurfaces[(i - _cw) * 8 + 0]
                            ->incrementTotalWt(wt, index + 3);
                    }
                    else
                    {
                        meshSurfaces[i * 8 + 0]
                            ->incrementTotalWt(wt2 / cosTheta, opposite_index);
                        meshSurfaces[i * 8 + 0]
                            ->incrementTotalWt(wt, opposite_index + 3);
                    }                  

                }
            }
        }
    }

    return;
}

void Solver::initializeWeights() 
{
    Track* track;
    int num_segments;
    std::vector<segment*> segments;
    segment* segment;
    int t, j, k, s;
    int num_threads = _num_azim / 2;
    MeshSurface **meshSurfaces = _geom->getMesh()->getSurfaces();

    /* Loop over each thread */
    for (t = 0; t < num_threads; t++) 
    {

        /* Loop over the pair of azimuthal angles for this thread */
        j = t;
        while (j < _num_azim) 
        {
            /* Loop over all tracks for this azimuthal angles */
            for (k = 0; k < _num_tracks[j]; k++) 
            {
                /* Initialize local pointers to important data structures */
                track = &_tracks[j][k];
                segments = track->getSegments();
                num_segments = track->getNumSegments();

                /* Loop over each segment in forward direction */
                for (s = 0; s < num_segments; s++) 
                {
                    segment = segments.at(s);
                    tallyLooWeight(track, segment, meshSurfaces, 1);
                }

                /* Loops over each segment in the reverse direction */
                for (s = num_segments-1; s > -1; s--) 
                {
                    segment = segments.at(s);
                    tallyLooWeight(track, segment, meshSurfaces, -1);
                }
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

    log_printf(NORMAL, "Number of corners tallied %d", _num_crn);
	
    MeshCell *meshCell;
    double w, h;
    for (int i = 0; i < _cw * _ch; i++)
    {
        meshCell = _geom->getMesh()->getCells(i);
        for (int s = 8 * i; s < 8 * i + 4; s++)
        {
            meshSurfaces[s]->setTotalWt((meshSurfaces[s]->getTotalWt(0) + 
                                         meshSurfaces[s]->getTotalWt(1)), 2);
            meshSurfaces[s]->setTotalWt((meshSurfaces[s]->getTotalWt(3) + 
                                         meshSurfaces[s]->getTotalWt(4)), 5);
            log_printf(NORMAL, "surface %d %.10f + %.10f = %.10f,"
                       " current-wt %.10f + %.10f = %.10f", s, 
                       meshSurfaces[s]->getTotalWt(0),
                       meshSurfaces[s]->getTotalWt(1),                   
                       meshSurfaces[s]->getTotalWt(2),
                       meshSurfaces[s]->getTotalWt(3),
                       meshSurfaces[s]->getTotalWt(4), 
                       meshSurfaces[s]->getTotalWt(5));
        }

        w = meshSurfaces[i * 8 + 1]->getTotalWt(2);
        h = meshSurfaces[i * 8 + 2]->getTotalWt(2);
        meshCell->setATWidth(w);
        meshCell->setATHeight(h);
        meshCell->setATL(0.5 * sqrt(w * w + h * h));
    }

    return;
} 

/* Performs MOC sweep(s), could be just one sweep or till convergance */
void Solver::MOCsweep(int max_iterations, int moc_iter) 
{
    Track* track;
    int num_segments;
    std::vector<segment*> segments;
    double* weights;
    segment* segment;
    double* polar_fluxes;
    double* scalar_flux;
    double* sigma_t;
    FlatSourceRegion* fsr;
    double fsr_flux[NUM_ENERGY_GROUPS];
    double* ratios;
    double delta;
    double volume;
    int t, j, k, s, p, e, pe;
    int num_threads = _num_azim / 2;
    MeshSurface **meshSurfaces = _geom->getMesh()->getSurfaces();

#if !STORE_PREFACTORS
    double sigma_t_l;
    int index;
#endif

    log_printf(INFO, "Fixed source iteration with max_iterations = %d and "
               "# threads = %d", max_iterations, num_threads);

    _num_crn = 0;

    int tally; 
    if ((_run_loo) && (!(_diffusion && (moc_iter == 0))))
    {
        tally = 2;
        log_printf(ACTIVE, "tally quadrature current");
    }
    else if ((_run_cmfd) || (_run_loo & _diffusion && (moc_iter == 0)))
    {
        tally = 1;
        log_printf(ACTIVE, "tally partial current");
    }
    else
        tally = 0;

    /* Loop for until converged or max_iterations is reached */
    for (int i = 0; i < max_iterations; i++)
    {
        /* Initialize flux in each region to zero */
        zeroFSRFluxes();
        zeroLeakage();

        /* Initializes mesh cells if ANY acceleration is on */
        if (_run_cmfd || _run_loo)
            zeroMeshCells();

        /* Loop over azimuthal each thread and azimuthal angle*
         * If we are using OpenMP then we create a separate thread
         * for each pair of reflecting azimuthal angles - angles which
         * wrap into cycles on each other */
#if USE_OPENMP && STORE_PREFACTORS
#pragma omp parallel for num_threads(num_threads)               \
    private(t, k, j, i, s, p, e, pe, track, segments,           \
            num_segments, weights, polar_fluxes,                \
            segment, fsr, ratios, delta, fsr_flux, currents)
#elif USE_OPENMP && !STORE_PREFACTORS
#pragma omp parallel for num_threads(num_threads)       \
    private(t, k, j, i, s, p, e, pe, track, segments,   \
            num_segments, weights, polar_fluxes,        \
            segment, fsr, ratios, delta, fsr_flux,      \
            sigma_t_l, index, currents)
#endif
        /* Loop over each thread */
        for (t = 0; t < num_threads; t++) 
        {

            /* Loop over the pair of azimuthal angles for this thread */
            j = t;
            while (j < _num_azim) 
            {
                /* Loop over all tracks for this azimuthal angles */
                for (k = 0; k < _num_tracks[j]; k++) 
                {
                    /* Initialize local pointers to important data structures */
                    track = &_tracks[j][k];
                    segments = track->getSegments();
                    num_segments = track->getNumSegments();
                    weights = track->getPolarWeights();
                    polar_fluxes = track->getPolarFluxes();

                    /* Loop over each segment in forward direction */
                    for (s = 0; s < num_segments; s++) 
                    {
                        segment = segments.at(s);
                        fsr = &_flat_source_regions[segment->_region_id];
                        ratios = fsr->getRatios();

                        /* Zero out temporary FSR flux array */
                        for (e = 0; e < NUM_ENERGY_GROUPS; e++)
                            fsr_flux[e] = 0.0;

                        /* Initialize the polar angle, energy group counter */
                        pe = 0;

#if !STORE_PREFACTORS
                        sigma_t = segment->_material->getSigmaT();

                        for (e = 0; e < NUM_ENERGY_GROUPS; e++)
                        {

                            fsr_flux[e] = 0;
                            sigma_t_l = sigma_t[e] * segment->_length;
                            sigma_t_l = std::min(sigma_t_l,10.0);
                            index = sigma_t_l / _pre_factor_spacing;
                            index = std::min(index * 2 * NUM_POLAR_ANGLES,
                                             _pre_factor_max_index);
							
                            for (p = 0; p < NUM_POLAR_ANGLES; p++)
                            {
                                delta = (polar_fluxes[pe] - ratios[e]) *
                                    (1 - (_pre_factor_array[index + 2 * p] 
                                          * sigma_t_l + 
                                          _pre_factor_array[index + 2 * p + 1]));
                                fsr_flux[e] += delta * weights[p];
                                polar_fluxes[pe] -= delta;
                                pe++;
								
                            }
                        }

#else
                        /* Loop over all polar angles and energy groups */
                        for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                        {
                            for (p = 0; p < NUM_POLAR_ANGLES; p++) 
                            {
                                delta = (polar_fluxes[pe] -ratios[e]) *
                                    segment->_prefactors[e][p];
                                fsr_flux[e] += delta * weights[p];
                                polar_fluxes[pe] -= delta;
                                pe++;
                            }
                        }
#endif

                        /* if segment crosses a surface in fwd direction, 
                           tally current/weight */
                        if (tally == 2)
                            tallyLooCurrent(track, segment, meshSurfaces, 1);
                        else if (tally == 1)
                            tallyCmfdCurrent(track, segment, meshSurfaces, 1);

                        /* Increments the scalar flux for this FSR */
                        fsr->incrementFlux(fsr_flux);
                    }

                    /* Transfers flux to outgoing track */
                    track->getTrackOut()->setPolarFluxes(track->isReflOut(), 0, 
                                                         polar_fluxes);

                    /* Tallys leakage */
                    pe = 0;
                    for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                    {
                        for (p = 0; p < NUM_POLAR_ANGLES; p++)
                        {
                            _geom->getSurface(track->getSurfFwd())
                                ->incrementLeakage
                                (track->isReflOut(), 
                                 polar_fluxes[pe]*weights[p] / 2.0, e);
                            pe++;
                        }
                    }

                    /* Loops over each segment in the reverse direction */
                    for (s = num_segments-1; s > -1; s--) 
                    {
                        segment = segments.at(s);
                        fsr = &_flat_source_regions[segment->_region_id];
                        ratios = fsr->getRatios();

                        /* Zero out temporary FSR flux array */
                        for (e = 0; e < NUM_ENERGY_GROUPS; e++)
                            fsr_flux[e] = 0.0;


                        /* Initialize the polar angle, energy group counter */
                        pe = GRP_TIMES_ANG;

#if !STORE_PREFACTORS
                        sigma_t = segment->_material->getSigmaT();
						
                        for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                        {
                            fsr_flux[e] = 0;
                            sigma_t_l = sigma_t[e] * segment->_length;
                            sigma_t_l = std::min(sigma_t_l,10.0);
                            index = sigma_t_l / _pre_factor_spacing;
                            index = std::min(index * 2 * NUM_POLAR_ANGLES,
                                             _pre_factor_max_index);

                            for (p = 0; p < NUM_POLAR_ANGLES; p++)
                            {
                                delta = (polar_fluxes[pe] - ratios[e]) *
                                    (1 - (_pre_factor_array[index + 2 * p] 
                                          * sigma_t_l
                                          + _pre_factor_array[index + 2 * p + 1]));
                                fsr_flux[e] += delta * weights[p];
                                polar_fluxes[pe] -= delta;
                                pe++;
                            }
                        }
#else
                        /* Loop over all polar angles and energy groups */
                        for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                        {
                            for (p = 0; p < NUM_POLAR_ANGLES; p++) 
                            {
                                delta = (polar_fluxes[pe] - ratios[e]) *
                                    segment->_prefactors[e][p];
                                fsr_flux[e] += delta * weights[p];
                                polar_fluxes[pe] -= delta;
                                pe++;
                            }
                        }
#endif

                        /* if segment crosses a surface in bwd direction, 
                           tally quadrature flux for LOO acceleration */
                        if (tally == 2)
                            tallyLooCurrent(track, segment, meshSurfaces, -1);
                        else if (tally == 1)
                            tallyCmfdCurrent(track, segment, meshSurfaces, -1);
						
                        /* Increments the scalar flux for this FSR */
                        fsr->incrementFlux(fsr_flux);
                    }

                    /* Transfers flux to incoming track */
                    track->getTrackIn()->setPolarFluxes(track->isReflIn(), 
                                                        GRP_TIMES_ANG, 
                                                        polar_fluxes);
                    /* Tallies leakage */
                    pe = GRP_TIMES_ANG;
                    for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                    {
                        for (p = 0; p < NUM_POLAR_ANGLES; p++)
                        {
                            _geom->getSurface(track->getSurfBwd())
                                ->incrementLeakage(track->isReflIn(), 
                                                   polar_fluxes[pe]*weights[p]
                                                   /2.0,
                                                   e);
                            pe++;
                        }
                    }
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

        /* If more than one iteration is requested, we only computes source for
         * the last iteration, all previous iterations are considered to be 
         * converging boundary fluxes */
        if (i == max_iterations - 1)
        {
#if USE_OPENMP
#pragma omp parallel for private(fsr, scalar_flux, ratios,  sigma_t, volume)
#endif
            /* Add in source term and normalize flux to volume for each FSR */
            for (int r = 0; r < _num_FSRs; r++) 
            {
                fsr = &_flat_source_regions[r];
                scalar_flux = fsr->getFlux();
                ratios = fsr->getRatios();
                sigma_t = fsr->getMaterial()->getSigmaT();
                volume = fsr->getVolume();

                for (int e = 0; e < NUM_ENERGY_GROUPS; e++) 
                {
                    fsr->setFlux(e, FOUR_PI * ratios[e] 
                                 + (scalar_flux[e] / 
                                    (2.0 * sigma_t[e] * volume)));
                }
            }
		
            normalizeFlux();

            /* computes new _k_eff; it is important that we compute new k 
             * before computing new source */
            _k_eff = computeKeff(i);

            /* Normalize scalar fluxes and computes Q for each FSR */
            updateSource();

            storeFsrFluxPower();

            /* Book-keeping: update old source in each FSR */
            double *source, *old_source;
            for (int r = 0; r < _num_FSRs; r++) 
            {
                fsr = &_flat_source_regions[r];
                source = fsr->getSource();
                old_source = fsr->getOldSource();
				
                for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
                    old_source[e] = source[e];
            }
        }
    } /* exit iteration loops */
		
    return;
} /* end of MOCsweep */

/* Normalizes the scalar flux in each FSR and angular flux track to that total 
 * volume-integrated fission source add up to the total volume. */
void Solver::normalizeFlux()
{			
    double fission_source = 0;
    double factor, volume;
    double* nu_sigma_f;
    double* scalar_flux;
    FlatSourceRegion* fsr;
    Material* material;
    int start_index, end_index;

    /* Compute total fission source for this region */
    for (int r = 0; r < _num_FSRs; r++) 
    {
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
    factor = _total_vol / fission_source;

#if USE_OPENMP
#pragma omp parallel for
#endif
    for (int r = 0; r < _num_FSRs; r++)
        _flat_source_regions[r].normalizeFluxes(factor);

    /* Renormalize angular boundary fluxes for each track */
#if USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++)
            _tracks[i][j].normalizeFluxes(factor);
    }

    /* Renormalize tallied current on each surface */
    if ((_run_cmfd) && !(_acc_after_MOC_converge))
    {
        int ng = NUM_ENERGY_GROUPS;
        MeshCell *meshCell;

        for (int i = 0; i < _cw * _ch; i++)
        {
            meshCell = _geom->getMesh()->getCells(i);
            for (int s = 0; s < 8; s++)
            {
                for (int e = 0; e < ng; e++)
                {
                    meshCell->getMeshSurfaces(s)->updateCurrent(factor, e);
                }
            }
        }
    }
    else if ((_run_loo) && !(_acc_after_MOC_converge))
    {
        int ng = NUM_ENERGY_GROUPS;
        MeshCell *meshCell;

        for (int i = 0; i < _cw * _ch; i++)
        {
            meshCell = _geom->getMesh()->getCells(i);
            for (int s = 0; s < 8; s++)
            {
                for (int e = 0; e < ng; e++)
                {
                    for (int ind = 0; ind < 2; ind++)
                    {
                        meshCell->getMeshSurfaces(s)->updateQuadCurrent(
                            factor, e, ind);
                    }
                }
            }
        }		
    }

    return;
}

/* Compute the source for each region */
void Solver::updateSource()
{
    double scatter_source, fission_source = 0;
    double* nu_sigma_f;
    double* sigma_s;
    double* chi;
    double* scalar_flux;
    double* source;
    FlatSourceRegion* fsr;
    Material* material;
    int start_index, end_index;

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
                scatter_source += sigma_s[G * NUM_ENERGY_GROUPS + g]
                    * scalar_flux[g];

            /* Set the total source for region r in group G */
            source[G] = ((1.0 / _k_eff) * fission_source *
                         chi[G] + scatter_source) * ONE_OVER_FOUR_PI;
        }
    }

    /* Update pre-computed source / sigma_t ratios */
    computeRatios();
}


double Solver::runLoo(int moc_iter)
{
    double loo_keff;

    _cmfd->storePreMOCMeshSource(_flat_source_regions);
    MOCsweep(_boundary_iteration + 1, moc_iter);
    _k_half = computeKeff(100);

    _cmfd->computeXS();
			 
    if (moc_iter == 0 && _diffusion == true)
    {
        _cmfd->computeDs();
        loo_keff = _cmfd->computeCMFDFluxPower(DIFFUSION, moc_iter);
    }
    else
    {
        _cmfd->computeQuadFlux();
        _cmfd->computeQuadSrc();
        loo_keff = _cmfd->computeLooFluxPower(moc_iter, _k_eff);
    }

    return loo_keff;
}

double Solver::runCmfd(int moc_iter)
{
    double cmfd_keff;

    MOCsweep(_boundary_iteration + 1, moc_iter);
    _k_half = computeKeff(100);

    /* compute cross sections and diffusion coefficients */
    //_cmfd->computeCurrent();
    _cmfd->computeXS();
    _cmfd->computeDs();

    /* Check for neutron balance */
    if (moc_iter == 10000)
        checkNeutronBalance();

    /* Run diffusion problem on initial geometry */
    if (moc_iter == 0 && _diffusion == true)
        cmfd_keff = _cmfd->computeCMFDFluxPower(DIFFUSION, moc_iter);
    else
        cmfd_keff = _cmfd->computeCMFDFluxPower(CMFD, moc_iter);

    return cmfd_keff;
}

double Solver::computeFsrL2Norm(double *old_fsr_powers)
{
    double l2_norm = 0.0;
    int num_counted = 0;
    for (int i = 0; i < _num_FSRs; i++)
    {
        if (old_fsr_powers[i] > 0.0)
        {
            log_printf(INFO, "new power = %f, old power = %f", 
                       _FSRs_to_powers[i], old_fsr_powers[i]);

            l2_norm += pow(_FSRs_to_powers[i] 
                           / old_fsr_powers[i] - 1.0, 2.0);
            num_counted++;
        }
    }
    l2_norm /= (double) num_counted;
    l2_norm = pow(l2_norm, 0.5);
   
    return l2_norm;
}

double Solver::computeFsrLinf(double *old_fsr_powers)
{
    double l2_norm = 0.0;
    double new_norm = 0.0;
    for (int i = 0; i < _num_FSRs; i++)
    {
        if (_FSRs_to_powers[i] > 0.0)
        {
            new_norm = fabs((double) old_fsr_powers[i] / 
                            ((double) _FSRs_to_powers[i]) - 1.0);
            l2_norm = fmax(new_norm, l2_norm);
        }
    }
    log_printf(DEBUG, "L2 norm = %e", l2_norm);
    return l2_norm;
}

double Solver::computeSpectralRadius(double *old_fsr_powers)
{
    double l2_norm = 0.0, rho = 0;
    int num_counted = 0; 
    //double tmp = 0, tmp_old = 0;
    for (int i = 0; i < _num_FSRs; i++)
    {
        if (old_fsr_powers[i] > 0.0)
        {
            log_printf(INFO, "new power = %f, old power = %f", 
                       _FSRs_to_powers[i], old_fsr_powers[i]);

            l2_norm += pow(_FSRs_to_powers[i] - old_fsr_powers[i], 2.0);
            num_counted++;
        }
    }
        /*
        tmp = 0;
        tmp_old = 0;
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            tmp += _FSRs_to_fluxes[e][i];
            tmp_old += old_fsr_fluxes[e][i];
        }
        l2_norm += pow(tmp - tmp_old, 2.0);
        */

    l2_norm /= (double) num_counted;
    l2_norm = pow(l2_norm, 0.5);

    _delta_phi.push(l2_norm);
    if (_delta_phi.size() == 3)
        _delta_phi.pop();

    rho = _delta_phi.back() / _delta_phi.front();
    
    return rho;
}

double Solver::kernel(int max_iterations) {
    int moc_iter;

    log_printf(NORMAL, "Starting kernel ...");

    /* Initial guess */
    _old_k_effs.push(_k_eff);
    _delta_phi.push(1.0);
    log_printf(NORMAL, "Starting guess of k_eff = %f", _k_eff);

    /* Check that each FSR has at least one segment crossing it */
    checkTrackSpacing();

    /* Check boundary conditions */
    checkBoundary();

    /* Set scalar flux to unity for each region */
    oneFSRFluxOldSource();
    initializeTrackFluxes(ONE_OVER_FOUR_PI);

    normalizeFlux();
    updateSource();

    double spectral_radius, eps_inf, eps_2;
    double old_fsr_powers[_num_FSRs];
    double **old_fsr_fluxes;
    old_fsr_fluxes = new double*[NUM_ENERGY_GROUPS];
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        old_fsr_fluxes[e] = new double[_num_FSRs];

    /* Computes and store initial FSR powers, store initial fluxes */
    storeFsrFluxPower();

    /* Stores the initial FSR powers into old_fsr_powers */
    for (int n = 0; n < _num_FSRs; n++)
    {
        old_fsr_powers[n] = _FSRs_to_powers[n];
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
            old_fsr_fluxes[e][n] = _FSRs_to_fluxes[e][n];
    }
    /* Source iteration loop */
    for (moc_iter = 0; moc_iter < max_iterations; moc_iter++) 
    {
        log_printf(INFO, "Iteration %d: k_eff = %f", moc_iter, _k_eff);

        /* Perform one sweep for no acceleration, or call one of the 
         * acceleration function which performs two sweeps plus acceleration */
        if ((_run_cmfd) && !(_acc_after_MOC_converge))
            _cmfd_k = runCmfd(moc_iter);
        else if ((_run_loo) && !(_acc_after_MOC_converge))
            _loo_k = runLoo(moc_iter);
        else 
            MOCsweep(1, moc_iter);

        /* Update FSR's flux based on cell-averaged flux coming from the
         * acceleration steps */
        if ((_run_cmfd || _run_loo) && !(_acc_after_MOC_converge))
            updateFlux(moc_iter);

        /* Computes the new keff */
        _k_eff = computeKeff(moc_iter);

        /* Computes new FSR source now we have new flux and new k */
        updateSource();

        /* We only store $k^{(m+1)}$; other intermediate keff does not matter */
        _old_k_effs.push(_k_eff);
        if (_old_k_effs.size() == NUM_KEFFS_TRACKED)
            _old_k_effs.pop();

        storeFsrFluxPower();

        spectral_radius = computeSpectralRadius(old_fsr_powers);
        eps_inf = computeFsrLinf(old_fsr_powers);
        eps_2 = computeFsrL2Norm(old_fsr_powers);

        /* book-keeping: save fsr fluxes and powers into old */
        for (int n = 0; n < _num_FSRs; n++)
        {
            old_fsr_powers[n] = _FSRs_to_powers[n];
            for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
                old_fsr_fluxes[e][n] = _FSRs_to_fluxes[e][n];
        }

        /* Prints out keff & eps, may update keff too based on _update_keff */
        printToScreen(moc_iter, eps_inf);
        printToLog(moc_iter, eps_inf, eps_2, spectral_radius);

        /* Alternative: if (_cmfd->getL2Norm() < _moc_conv_thresh) */
        if (eps_2 < _moc_conv_thresh) 
        {
            /* Run one steps of acceleration if it is requested to do so */
            if (_acc_after_MOC_converge)
            {
                if (_run_loo)
                    _loo_k = runLoo(10000);
                if (_run_cmfd)
                    _cmfd_k = runCmfd(10000);
            }

            printToMinimumLog(moc_iter);
            plotEverything(moc_iter);

            return _k_eff;
        }
    }

    log_printf(WARNING, "Unable to converge the source after %d iterations",
               max_iterations);

    return _k_eff;
}

void Solver::plotEverything(int moc_iter)
{
    /* plot pin powers */
    if (_compute_powers)
        plotPinPowers();

    /* plot CMFD flux and xs */
    if (_run_cmfd && _plotter->plotCurrent() )
    {
        _cmfd->computeXS();
        _cmfd->computeDs();
        _plotter->plotDHats(_geom->getMesh(), moc_iter);
        //_plotter->plotNetCurrents(_geom->getMesh(), moc_iter);
        _plotter->plotXS(_geom->getMesh(), moc_iter);
    }

    /* plot LOO flux and xs */
    if (_run_loo)
    {
        _plotter->plotQuadFlux(_geom->getMesh(), moc_iter);
        _plotter->plotNetCurrents(_geom->getMesh(), moc_iter);
        _plotter->plotXS(_geom->getMesh(), moc_iter);
    }

    /* plot FSR scalar flux */
    if (_plotter->plotFlux())
        plotFluxes(moc_iter);

    return;
}

void Solver::printToMinimumLog(int moc_iter)
{
    std::ofstream logfile;
    std::stringstream string;
    string << "benchmark.txt";
    std::string title_str = string.str();
    logfile.open(title_str.c_str(), std::fstream::app);
    logfile << _geometry_file << " " 
            <<  std::setprecision(11) <<  _k_eff;
    logfile << " " << moc_iter << std::endl;
    logfile.close();

    return;
}


void Solver::printToLog(int moc_iter, double eps_inf, double eps_2, double rho)
{
    std::ofstream logfile;
    std::stringstream string;

    string << "l2_norm_" << (_num_azim*2) 
           << "_" << std::fixed 
           << std::setprecision(2) <<  _track_spacing
           << "_bi_" << _boundary_iteration
           << "_"
           << std::setprecision(1) << _damp_factor;

    if (_update_boundary)
        string << "_update";
    else
        string << "_noupda";

    if (_run_cmfd)
        string << "_cmfd.txt";
    else if (_run_loo1)
        string << "_loo1.txt";
    else if (_run_loo2)
        string << "_loo2.txt";
    else
        string << "_unac.txt";


    std::string title_str = string.str();

    if (moc_iter == 0)
    {
        logfile.open(title_str.c_str(), std::fstream::trunc);
        logfile << "# iteration,"
                << " cell l2 norm (m+1/2, m+1),"
                << " fsr l-inf norm (m, m+1),"
                << " fsr l2 norm (m, m+1)," 
                << " keff relative change,"
                << " #lo iterations, "
                << " keff,"
                << " spectral radius"
                << std::endl;
    }
    else
    {
        logfile.open(title_str.c_str(), std::ios::app);
        logfile << moc_iter 
                << " " << _cmfd->getL2Norm() 
                << " " << eps_inf
                << " " << eps_2
                << " " << 1.0 -  _old_k_effs.front() / _old_k_effs.back()
                << " " << _cmfd->getNumIterToConv() 
                << " " << std::setprecision(11) << _old_k_effs.back()
                << " " << rho
                << std::endl;
    }

    logfile.close();
        
}



/*
 * check neutron balance in each mesh cell
 */
void Solver::checkNeutronBalance()
{
    log_printf(INFO, "Checking neutron balance...");
    double leak = 0, absorb = 0, fis = 0, src = 0;
    double tot_leak = 0, tot_absorb = 0, tot_fis = 0, tot_src = 0;
    double flux, vol, residual;
    MeshCell *meshCell, *meshCellNext;
    MeshSurface* surf;
    Mesh *mesh = _geom->getMesh();
    int cell_height = mesh->getCellHeight();
    int cell_width = mesh->getCellWidth();
    int ng = NUM_ENERGY_GROUPS;
    if (mesh->getMultigroup() == false)
        ng = 1;

    /* loop over mesh cells in y direction */
    for (int y = 0; y < cell_height; y++)
    {
        for (int x = 0; x < cell_width; x++)
        {
            leak = 0; 
            absorb = 0;
            fis = 0;
            src = 0;

            /* get mesh cell */
            meshCell = mesh->getCells(y * cell_width + x);

            for (int e = 0; e < ng; e++)
            {
                leak = 0;
                absorb = 0;
                fis = 0;
                src = 0;

                flux = meshCell->getOldFlux()[e];
                vol = meshCell->getATVolume();
				
                absorb += meshCell->getSigmaA()[e] * flux;

                double tmp_fis = 0.0;
                for (int g = 0; g < ng; g++)
                {
                    tmp_fis += meshCell->getNuSigmaF()[g] * 
                        meshCell->getOldFlux()[g];
                }
                fis += meshCell->getChi()[e] * tmp_fis;

                for (int g = 0; g < ng; g++)
                {
                    if (g != e)
                    {
                        src += meshCell->getSigmaS()[g * ng + e] * 
                            meshCell->getOldFlux()[g];
                    }
                }
                src += fis / _k_eff;

                for (int s = 0; s < 4; s++)
                {
                    surf = meshCell->getMeshSurfaces(s);
                    leak += surf->getCurrent(e);
                }

                if (x > 0)
                {
                    meshCellNext = mesh->getCells(y * cell_width + x - 1);
                    leak -= meshCellNext->getMeshSurfaces(2)->getCurrent(e);
                }
                else
                    leak -= meshCell->getMeshSurfaces(0)->getCurrent(e);
				
                if (x < cell_width - 1)
                {
                    meshCellNext = mesh->getCells(y * cell_width + x + 1);
                    leak -= meshCellNext->getMeshSurfaces(0)->getCurrent(e);
                }
                else
                    leak -= meshCell->getMeshSurfaces(2)->getCurrent(e);

                if (y > 0)
                {
                    meshCellNext = mesh->getCells((y - 1) * cell_width + x);
                    leak -= meshCellNext->getMeshSurfaces(1)->getCurrent(e);
                }
                else
                    leak -= meshCell->getMeshSurfaces(3)->getCurrent(e);

                if (y < cell_height - 1)
                {
                    meshCellNext = mesh->getCells((y + 1) * cell_width + x);
                    leak -= meshCellNext->getMeshSurfaces(3)->getCurrent(e);
                }
                else
                    leak -= meshCell->getMeshSurfaces(1)->getCurrent(e);

                leak /= vol;

                /* compute total residual and average ratio */
                /* residual = leakage + absorption - fission */
                residual = leak + absorb - src;
                log_printf(ACTIVE, "CMFD cell %d energy %d, residual %.10f"
                           " leak: %.10f"
                           " absorb: %.10f, src: %.10f,"
                           "  fis: %.10f, keff: %.10f", 
                           y * cell_width + x, e, 
                           residual,
                           leak,
                           absorb, src, 
                           fis, fis / (leak + absorb));
                tot_leak += leak;
                tot_absorb += absorb;
                tot_fis += fis;
                tot_src += src;
            }	
        }
    }

    log_printf(NORMAL, "CMFD over all cell, keff: %.10f"
               " fis: %f, absorb: %f, leak: %f, src = %f, res = %f", 
               tot_fis / (tot_leak + tot_absorb),
               tot_fis, tot_absorb, tot_leak, 
               tot_src, tot_leak + tot_absorb - tot_src);
	
    return;
}



void Solver::checkNeutronBalanceWithDs()
{
    log_printf(INFO, "Checking neutron balance...");
    MeshCell *meshCell;
    Mesh *mesh = _geom->getMesh();
    double leak = 0, absorb = 0, fis = 0;
    int cell_height = mesh->getCellHeight();
    int cell_width = mesh->getCellWidth();
    int ng = NUM_ENERGY_GROUPS;
    if (mesh->getMultigroup() == false)
        ng = 1;

    /* loop over mesh cells in y direction */
    for (int xy = 0; xy < cell_height * cell_width; xy++)
    {
        leak = 0; 
        absorb = 0;
        fis = 0;

        /* get mesh cell */
        meshCell = mesh->getCells(xy);

        /* leakage */
        for (int s = 0; s < 4; s++)
        {
            for (int e = 0; e < ng; e++)
            {
                if (s == 0)
                {
                    leak += meshCell->getMeshSurfaces(s)->getDHat()[e] 
                        * meshCell->getOldFlux()[e]
                        * meshCell->getHeight();
                    leak += meshCell->getMeshSurfaces(s)->getDTilde()[e]
                        * meshCell->getOldFlux()[e]
                        * meshCell->getHeight();
                }
                else if (s == 2){
                    leak += meshCell->getMeshSurfaces(s)->getDHat()[e] 
                        * meshCell->getOldFlux()[e]
                        * meshCell->getHeight();
                    leak -= meshCell->getMeshSurfaces(s)->getDTilde()[e]
                        * meshCell->getOldFlux()[e]
                        * meshCell->getHeight();
                }
                else if (s == 1){
                    leak += meshCell->getMeshSurfaces(s)->getDHat()[e] 
                        * meshCell->getOldFlux()[e]
                        * meshCell->getWidth();
                    leak -= meshCell->getMeshSurfaces(s)->getDTilde()[e]
                        * meshCell->getOldFlux()[e]
                        * meshCell->getWidth();
                }
                else if (s == 3){
                    leak += meshCell->getMeshSurfaces(s)->getDHat()[e] 
                        * meshCell->getOldFlux()[e]
                        * meshCell->getWidth();
                    leak += meshCell->getMeshSurfaces(s)->getDTilde()[e]
                        * meshCell->getOldFlux()[e]
                        * meshCell->getWidth();
                }
            }
        }
			
        /* compute absorb and fis rates */
        for (int e = 0; e < ng; e ++)
        {
            absorb += meshCell->getSigmaA()[e] * meshCell->getATVolume() * 
                meshCell->getOldFlux()[e];
            fis += meshCell->getNuSigmaF()[e] * meshCell->getATVolume()
                * meshCell->getOldFlux()[e];
        }	

        /* compute total residual and average ratio */
        /* residual = leakage + absorption - fission */
        double residual = leak + absorb - fis;
        log_printf(NORMAL, "CMFD residual %.10f"
                   " fis: %.10f, absorb: %.10f, leak: %.10f, keff: %.10f", 
                   residual, fis, absorb, leak, fis / (leak + absorb));
    } /* end of looping through cells */

    return;
}


FlatSourceRegion* Solver::getFSRs(){

    return _flat_source_regions;
}


/* Set the Old Flux for each FSR equal to FSR flux */
void Solver::setOldFSRFlux()
{
    FlatSourceRegion* fsr;

    /* Compute total fission source for this region */
    for (int r = 0; r < _num_FSRs; r++) 
    {
        /* Get pointers to important data structures */
        fsr = &_flat_source_regions[r];

        /* loop over energy groups */
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            fsr->setOldFlux(e, fsr->getFlux()[e]);

        }
    }
}

