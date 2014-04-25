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

    _first_diffusion = opts->getFirstDiffusion();
    _num_first_diffusion = opts->getNumFirstDiffusion();
    _acc_after_MOC_converge = opts->getAccAfterMOCConverge();
    _k_eff = opts->getKGuess();
    _damp_factor = opts->getDampFactor();
    _track_spacing = opts->getTrackSpacing();
    _boundary_iteration = opts->getBoundaryIteration();
    _update_boundary = opts->getUpdateBoundary();
    _use_up_scattering_xs = opts->getUseUpScatteringXS();

    _geometry_file = opts->getGeometryFile();
    _geometry_file_no_slash = opts->getGeometryFile();
    for (unsigned i = 0; i < _geometry_file.length(); i++)
    {
        switch(_geometry_file[i]) {
        case '/':
            _geometry_file_no_slash[i] = '_';
        }
    }
    std::stringstream string;

    string << _geometry_file_no_slash << "_" << (_num_azim*2) 
           << "_" << std::fixed 
           << std::setprecision(2) <<  _track_spacing
        //<< "_bi_" << _boundary_iteration
           << "_"
           << std::setprecision(1) << _damp_factor;

    if (_update_boundary)
        string << "_update";
    else
        string << "_noupda";

    if (_use_up_scattering_xs)
        string << "_upscat";
    else
        string << "_noupsc";

    if (_run_cmfd)
        string << "_cmfd.txt";
    else if (_run_loo1)
        string << "_loo1.txt";
    else if (_run_loo2)
        string << "_loo2.txt";
    else
        string << "_unac.txt";

    _log_file = string.str();

    _reflect_outgoing = opts->getReflectOutgoing();

    _nq = 4;
    
    _plot_loo = opts->plotQuadFlux();
    _plot_flux = opts->plotFluxes();

    _acc_k = _k_eff;

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

    _cw = _geom->getMesh()->getCellWidth();
    _ch = _geom->getMesh()->getCellHeight();

    _plotter->setFSRs(_flat_source_regions);

    /* Gives cmfd a pointer to the FSRs */
    if (_run_cmfd || _run_loo || _acc_after_MOC_converge) 
    {
        _cmfd->setFSRs(_flat_source_regions);
        _cmfd->setTracks(_tracks);

        initializeWeights();
        
        log_printf(NORMAL, "Acceleration on:"
                   " upscattering = %d, boundary iteration = %d,"
                   " damping = %.2f, update boundary flux = %d,"
                   " reflect outgoing = %d", 
                   _use_up_scattering_xs, _boundary_iteration, 
                   _damp_factor, _update_boundary, _reflect_outgoing);
    }

    _pin_powers.assign(_ch * _cw, 1.0);
}

void Solver::runCmfd() {
    _run_cmfd = true;
    _run_loo = false;
    _run_loo1 = false;
    _run_loo2 = false;
    _cmfd->runCmfd();
}

void Solver::runLoo1() {
    _run_cmfd = false;
    _run_loo = true;
    _run_loo1 = true;
    _run_loo2 = false;
    _cmfd->runLoo1();
}

void Solver::runLoo2() {
    _run_cmfd = false;
    _run_loo = true;
    _run_loo1 = false;
    _run_loo2 = true;
    _cmfd->runLoo2();
}

void Solver::setK(double k) {
    _k_eff = k;
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
    log_printf(INFO, "Initializing FSRs...");

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
    for (int r = 0; r < _num_FSRs; r++) {
        /* Set the id */
        _flat_source_regions[r].setId(r);

        /* Get the cell corresponding to this FSR from the geometry */
        cell = static_cast<CellBasic*>(_geom->findCell(univ_zero, r));

        /* Get the cell's material and assign it to the FSR */
        material = _geom->getMaterial(cell->getMaterial());
        _flat_source_regions[r].setMaterial(material);
        _flat_source_regions[r].setQuadId(cell->getQuadId());

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

    log_printf(INFO, "Setting all non-vacuum boundary fluxes to %f...", flux);

    double* polar_fluxes;

    log_printf(DEBUG, "geometry: %f %f", _geom->getWidth(),  
               _geom->getWidth());

    /* Loop over azimuthal angle, track, polar angle, energy group
     * and set each track's incoming and outgoing flux to zero */
    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++) 
        {
            polar_fluxes = _tracks[i][j].getPolarFluxes();
            /* Be careful: (x,y) spans from - width to + width */
            double x0 = _tracks[i][j].getStart()->getX();
            double y0 = _tracks[i][j].getStart()->getY();
            double x1 = _tracks[i][j].getEnd()->getX();
            double y1 = _tracks[i][j].getEnd()->getY();

            /* forward */
            if (onVacuumBoundary(x0, y0))
            {
                log_printf(DEBUG, "found vacuum at (%f %f)", x0, y0);
                for (int i = 0; i < GRP_TIMES_ANG; i++)
                    polar_fluxes[i] = 0.0;
            }
            else
            {
                for (int i = 0; i < GRP_TIMES_ANG; i++)
                    polar_fluxes[i] = flux;
            }

            /* backward */
            if (onVacuumBoundary(x1, y1))
            {
                log_printf(DEBUG, "found vacuum at (%f %f)", x1, y1);
                for (int i = GRP_TIMES_ANG; i < 2 * GRP_TIMES_ANG; i++)
                    polar_fluxes[i] = 0.0;
            }
            else
            {
                for (int i = GRP_TIMES_ANG; i < 2 * GRP_TIMES_ANG; i++)
                    polar_fluxes[i] = flux;
            }
        }
    }
    return;
}

bool Solver::onVacuumBoundary(double x, double y)
{
    for (int s = 0; s < 4; s++)
    {
        if ((_geom->getMesh()->getBoundary(s) == VACUUM) && onBoundary(x, y, s))
            return true;
    }
    return false;
}

bool Solver::onBoundary(double x, double y, int s)
{
    double w = _geom->getWidth() / 2.0;
    double l = _geom->getHeight() / 2.0;
    double delta = 1e-5;

    if ((s == 0) && (x < -w + delta))
        return true;
    else if ((s == 1) && (y < -l + delta))
        return true;
    else if ((s == 2) && (x > w - delta))
        return true;
    else if ((s == 3) && (y > l - delta))
        return true;

    return false;
}


/**
 * Set the scalar flux for each energy group inside each FSR
 * to unity
 */
void Solver::oneFSRFlux() 
{
    log_printf(INFO, "Setting all FSR scalar fluxes to unity...");
    FlatSourceRegion* fsr;

    /* Loop over all FSRs and energy groups */
    for (int r = 0; r < _num_FSRs; r++) 
    {
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

    log_printf(INFO, "Setting all mesh cell fluxes and currents to zero...");

    /* get mesh */
    Mesh* mesh = _geom->getMesh();
    MeshCell* meshCell;

    /* loop over mesh cells */
    for (int i = 0; i < _cw * _ch; i++)
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
                for (int j = 0; j < _nq; j++)
                    meshCell->getMeshSurfaces(surface)->setQuadCurrent(0, e, j);
            }
        }
    }
    return;
}


void Solver::zeroLeakage()
{
    for (int s = 1; s < 5; s++)
    {
        log_printf(DEBUG, "geometry surface %d has leakage %f", s, 
                   _geom->getSurface(s)->getLeakage()[0]);
    }

    for (int s = 1; s < 5; s++)
        _geom->getSurface(s)->zeroLeakage();
    return;
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

// _geom->getSurface() uses index 1 through 4. These leakage terms are
// surface integrated current. 
leakage = 0;
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

log_printf(DEBUG, "MOC k = %f / (%f + %f) = %f", 
           tot_fission, tot_abs, leakage, k);
//_geom->getMesh()->setKeffMOC(_k_eff, iteration);
return k;
}


/**
 * Update FSR's scalar fluxes, normalize them and update the source
 */
void Solver::prolongation(int moc_iter) 
{
    if (moc_iter > -10) //47 (2x2_leakage) 79 (2x2) 229 (4x4_leakage) 
    {
        log_printf(DEBUG, " iter %d scalar flux prolongation", moc_iter);
        _cmfd->updateFSRScalarFlux(moc_iter);

        if (_update_boundary)
        {
            /* standard LOO update */
            if (_run_loo && 
                (!(_first_diffusion && (moc_iter < _num_first_diffusion))))
            {
                log_printf(DEBUG, " iter %d boundary angular flux "
                           "prolongation: by quadrature",
                           moc_iter);
                updateBoundaryFluxByQuadrature(moc_iter);
                //_cmfd->updateBoundaryFluxByScalarFlux(moc_iter);
                //_cmfd->updateBoundaryFluxBySrc(moc_iter);
                //_cmfd->updateBoundaryFluxByNetCurrent(moc_iter);
            }
            /* first diffusion step or standard CMFD */
            else 
            {
                log_printf(DEBUG, " iter %d boundary angular flux "
                           "prolongation: by scalar flux",
                           moc_iter);
                _cmfd->updateBoundaryFluxByScalarFlux(moc_iter);
            }
        }
    }

    /* Debug routine to force all boundaries to go to zero */
    //zeroVacuumBoundaries();

    return;
}

void Solver::storeMOCBoundaryFlux()
{
    MeshSurface *meshSurface, **meshSurfaces;
    meshSurfaces = _geom->getMesh()->getSurfaces();
    for (int i = 0; i < 8 * _cw * _ch; i++)
    {
        meshSurface = meshSurfaces[i];
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            for (int ind = 0; ind < _nq; ind++)
            {
                meshSurface->setOldQuadFlux(meshSurface->getQuadFlux(e, ind),
                                        e, ind);
            }
        }
    }
}


void Solver::updateBoundaryFluxByQuadrature(int moc_iter)
{
    Track *track;
    segment *seg;
    MeshSurface *meshSurface1, *meshSurface2;
    MeshSurface **meshSurfaces;  
    double phi, factor;
    int ind, num_segments, pe, surf_id, num_updated = 0;
    int this_ind;
    meshSurfaces = _geom->getMesh()->getSurfaces();

    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++)
        {
            track = &_tracks[i][j];
            num_segments = track->getNumSegments();
            phi = track->getPhi();

            /* Use the opposite quadrature, say index 1 for < pi/2.
             * regardless of whether _reflect_outgoing is turned on or
             * not, this is always the case because during LOO perfect
             * reflectiveness is assumed. */
             if (phi < PI / 2.0) 
             {
                 ind = 1; 
                 this_ind = 0;
             }
             else 
             {
                 ind = 0;
                 this_ind = 0;
             }

            /* Forward direction is 0, backward is 1 */
            for (int dir = 0; dir < 2; dir++)
            {
                seg = track->getSegment(dir * (num_segments - 1));
         
                /* find the boundary surface: if forward direction,
                 * find the surface that the segment starts from.  */
                if (dir == 0)
                    surf_id = seg->_mesh_surface_bwd; 
                else
                    surf_id = seg->_mesh_surface_fwd;

                int corner_id = surf_id % 8;
                int cell_id = surf_id / 8;
                int y = cell_id / _cw;
                int x = cell_id % _cw; 

#if 0
                /* FIXME: need to comment out to get c5g7_cc working */
                if (moc_iter > -1)
                {
                if ((x == 0) && (_geom->getMesh()->getBoundary(0) == VACUUM) &&
                    (corner_id == 0))
                {
                    if (meshSurfaces[surf_id]->getOldQuadFlux(0,ind) > 1e-8)
                    {
                        log_printf(NORMAL, "cell %d (%d %d) surface %d has old"
                                   " flux %f new flux %f", 
                                   cell_id, x, y, corner_id, 
                                   meshSurfaces[surf_id]->getOldQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,ind));
                    }
                    break;
                }
                if ((x == _cw - 1) 
                    && (_geom->getMesh()->getBoundary(2) == VACUUM) &&
                    (corner_id == 2))
                    {
                        log_printf(NORMAL, "cell %d (%d %d) surface %d has old"
                                   " flux %f new flux %f, this index %f %f", 
                                   cell_id, x, y, corner_id, 
                                   meshSurfaces[surf_id]->getOldQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getOldQuadFlux(0, this_ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,this_ind)
                            );
                    }
                    break;
                if ((y == 0) && (_geom->getMesh()->getBoundary(3) == VACUUM) &&
                    (corner_id == 3))
                    {
                        log_printf(NORMAL, "cell %d (%d %d) surface %d has old"
                                   " flux %f new flux %f, this index %f %f", 
                                   cell_id, x, y, corner_id, 
                                   meshSurfaces[surf_id]->getOldQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getOldQuadFlux(0, this_ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,this_ind)
                            );
                    }
                    break;
                if ((y == _ch - 1) 
                    && (_geom->getMesh()->getBoundary(1) == VACUUM)
                    && (corner_id ==1))
                    {
                        log_printf(NORMAL, "cell %d (%d %d) surface %d has old"
                                   " flux %f new flux %f, this index %f %f", 
                                   cell_id, x, y, corner_id, 
                                   meshSurfaces[surf_id]->getOldQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,ind),
                                   meshSurfaces[surf_id]->getOldQuadFlux(0, this_ind),
                                   meshSurfaces[surf_id]->getQuadFlux(0,this_ind)
                            );
                    }
                    break;
                }
#endif

                int surf1 = surf_id, surf2 = surf_id;

                if (corner_id < 4)
                {
                    surf1 = surf_id;
                    surf2 = surf_id;
                }
                else if (corner_id < 5)
                {
                    if (x == 0) /* left column */
                    {
                        surf1 = surf_id - 4;
                        surf2 = surf1 + 8 * _cw;
                    }
                    else if (y == _ch - 1) /* bottom row */
                    {
                        surf1 = surf_id - 3;
                        surf2 = surf1 - 8;
                    }
                    else
                    {
                        log_printf(ERROR, "boundary update fail in cell %d"
                                   " (x: %d y: %d) corner %d", cell_id, 
                                   x, y, corner_id);
                    }
                }
                else if (corner_id < 6)
                {
                    if (x == _cw - 1) /* right column */
                    {
                        surf1 = surf_id - 3;
                        surf2 = surf1 + 8 * _cw;
                    }
                    else if (y == _ch - 1) /* bottom row */
                    {
                        surf1 = surf_id - 4;
                        surf2 = surf1 + 8;
                    }
                    else
                    {
                        log_printf(ERROR, "boundary update fail in cell %d"
                                   " (x: %d y: %d) corner %d", cell_id, 
                                   x, y, corner_id);
                    }
                }
                else if (corner_id < 7)
                {
                    if (x == _cw - 1) /* right column */
                    {
                        surf1 = surf_id - 4;
                        surf2 = surf1 - 8 * _cw;
                    }
                    else if (y == 0) /* top row */
                    {
                        surf1 = surf_id - 3;
                        surf2 = surf1 + 8;
                    }
                    else
                    {
                        log_printf(ERROR, "boundary update fail in cell %d"
                                   " (x: %d y: %d) corner %d", cell_id, 
                                   x, y, corner_id);
                    }
                }
                else
                {
                    if (x == 0) /* left column */
                    {
                        surf1 = surf_id - 7;
                        surf2 = surf1 - 8 * _cw;
                    }
                    else if (y == 0) /* top row */
                    {
                        surf1 = surf_id - 4;
                        surf2 = surf1 - 8;
                    }
                    else
                    {
                        log_printf(ERROR, "boundary update fail in cell %d"
                                   " (x: %d y: %d) corner %d", cell_id, 
                                   x, y, corner_id);
                    }
                }


                assert((surf1 % 8) < 4);
                assert(surf1 >= 0);
                assert(surf1 < 8 * _cw * _ch);
                assert((surf2 % 8) < 4);
                assert(surf2 >= 0);
                assert(surf2 < 8 * _cw * _ch);


                meshSurface1 = meshSurfaces[surf1];
                meshSurface2 = meshSurfaces[surf2];

                /* initialize starting point of pe which is 0 for
                 * forward direction and GRP_TIMES_ANG for backward
                 * direction */
                pe = dir * GRP_TIMES_ANG;

                for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
                {
                    if (meshSurface2->getOldQuadFlux(e, ind) < 1e-8)
                    {		
                        pe += NUM_POLAR_ANGLES;
                    }
                    /*
                    else if (meshSurface1->getOldQuadFlux(e,ind) < 1e-8)
                    {
                        factor = 0.5 * (meshSurface2->getQuadFlux(e, ind) + 
                                        meshSurface2->getQuadFlux(e, ind));
                        for (int p = 0; p < NUM_POLAR_ANGLES; p++)
                        {
                            track->setPolarFluxesByIndex(pe, factor);
                            pe++;
                            num_updated++;
                        }
                    }
                    */
                    else
                    {
                        factor = 0.5 * (meshSurface1->getQuadFlux(e, ind)
                                        / meshSurface1->getOldQuadFlux(e, ind)
                                        + meshSurface2->getQuadFlux(e, ind)
                                        / meshSurface2->getOldQuadFlux(e, ind));

                        factor = _damp_factor * factor + 1.0 - _damp_factor; 
                        for (int p = 0; p < NUM_POLAR_ANGLES; p++)
                        {
                            track->updatePolarFluxes(pe, factor);
                            pe++;
                            num_updated++;
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

/* Prints & Update keff for MOC sweep. Causion: prolongated results
 * are printed as moc_iter + 1, as we want to see on the screen iter
 * 1,... instead of 0-indexed results */
void Solver::printToScreen(int moc_iter)
{
    if ((moc_iter < _num_first_diffusion + 1) && _run_cmfd && _first_diffusion)
    {
        printf("Iter %d, MOC k^(m+1) = %.10f, Diffusion k = %.10f,"
               " MOC k^(m+1/2) = %.10f" 
               " FS eps = %.4e,  k eps = %.4e, #Diffusion = %d\n", 
               moc_iter, _k_eff, _acc_k, _k_half, _old_eps_2.back(),
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());
    }
    else if ((moc_iter < _num_first_diffusion + 1) && _run_loo 
             && _first_diffusion)
    {
        printf("Iter %d, MOC k^(m+1) = %.10f, Diffusion k = %.10f,"
               " MOC k^(m+1/2) = %.10f" 
               " FS eps = %.4e,  k eps = %.4e, #Diffusion = %d\n", 
               moc_iter, _k_eff, _acc_k, _k_half, _old_eps_2.back(),
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());
    }
    else if ((_run_cmfd) && !(_acc_after_MOC_converge))
    {
        printf("Iter %d, MOC k^(m+1) = %.10f,"
               " MOC k^(m+1/2) = %.10f" 
               " FS eps = %.2e, FS ratio = %.2f, k eps = %.2e, #CMFD = %d\n", 
               moc_iter, _k_eff, _k_half, _old_eps_2.back(),
               _old_eps_2.front() / _old_eps_2.back(), 
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());
    }
    else if ((_run_loo) && !(_acc_after_MOC_converge))
    {
        printf("Iter %d, MOC k^(m+1) = %.10f, LOO k = %.10f,"
               " MOC k^(m+1/2) = %.10f"
               " FS eps = %.4e, FS ratio = %f, k eps = %.4e, #LOO = %d\n", 
               moc_iter, _k_eff, _acc_k, _k_half,  _old_eps_2.back(),
               _old_eps_2.front() / _old_eps_2.back(), 
               (_old_k_effs.back() - _old_k_effs.front()) / _old_k_effs.back(),
               _cmfd->getNumIterToConv());
    }
    else
    {
        printf("Iter %d, MOC k = %.10f, FS eps = %.4e, "
               "FS ratio = %f, k eps = %.4e\n", 
               moc_iter, _k_eff,  _old_eps_2.back(),
               _old_eps_2.front() / _old_eps_2.back(), 
               (_old_k_effs.front() - _old_k_effs.back()) 
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

/* Print out the four boundary conditions around the exterior geometry in
 * the order of left, bottom, right, top.  
 */
void Solver::checkBoundary()
{
    std::string bc[4];
    for (int s = 0; s < 4; s++)
    {
        bc[s] = "unset";

        switch (_geom->getMesh()->getBoundary(s))
        {
        case REFLECTIVE:
            bc[s] = "reflective";
            break;
        case VACUUM:
            bc[s] = "vacuum";
            break;
        case BOUNDARY_NONE:
            bc[s] = "unknown";
            break;
        }
    }
    log_printf(NORMAL, "Boundary conditions = %s %s %s %s", 
               bc[0].c_str(), bc[1].c_str(), bc[2].c_str(), bc[3].c_str()); 
    return;
}
/**
 * Plot the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates
 */
void Solver::computePinPowers() {

    log_printf(INFO, "Computing pin powers...");

    FlatSourceRegion* fsr;

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

    for (int i = 0; i < _cw * _ch; i++)
    {
        int j = _geom->getMesh()->getCells(i)->getFSRStart();
        _pin_powers.push_front(_FSRs_to_pin_powers[j]);
    }

/*
    avg_pin_power = tot_pin_power / num_nonzero_pins;

    for (int i=0; i < _num_FSRs; i++) {
        _FSRs_to_pin_powers[i] /= avg_pin_power;
    }
*/
    return;
}

double Solver::computePinPowerNorm()
{
    double new_power, old_power, norm = 0; 
    int counter = 0;
    std::forward_list<double> old_pin_powers;

    /* copy all items from _pin_powers into old_pin_powers; clear the former */
    while (!_pin_powers.empty())
    {
        old_pin_powers.push_front(_pin_powers.front());
        _pin_powers.pop_front();
    }

    /* fill up _pin_powers with new values */
    computePinPowers();
    
    for (auto it = _pin_powers.begin(); it != _pin_powers.end(); it++)
    {
        new_power = *it;
        old_power = old_pin_powers.front();
        if (new_power > 1e-10)
        {
            log_printf(DEBUG, "new power = %e",
                       new_power - 2.25);
            norm += pow(new_power / old_power - 1.0, 2);
            counter += 1;
        }
        old_pin_powers.pop_front();
    }

    norm /= (double) counter;
    norm = sqrt(norm);
    return norm;
}

/**
 * Plot the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates
 */
void Solver::plotPinPowers() {
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

    computePinPowers();

    /* make FSR BitMap */
    _plotter->copyFSRMap(bitMapFSR->pixels);

    log_printf(NORMAL, "Plotting pin powers...");
    _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                            _FSRs_to_pin_powers);
    plot(bitMap, "pin_powers", _plotter->getExtension());

    /* delete bitMaps */
    deleteBitMap(bitMapFSR);
    deleteBitMap(bitMap);

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
    for (int i = 0; i < 1; i++)
    {
        string.str("");
        string << "iter" << moc_iter << "flux" << i + 1 << "group";
        title_str = string.str();
        _plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels, 
                                //_FSRs_to_fluxes[i]);
                                _cmfd->getCellSource());
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


void Solver::tallyLooCurrentIncoming(Track *track, segment *segment, 
                                     MeshSurface **meshSurfaces, int direction)
{
    int index = 2, pe = 0, pe_initial = 0, p, e, surfID;
    MeshSurface *meshSurface;
    double *weights, *polar_fluxes;

    /* Get the ID of the surface that the segment starts on (forward),
     * ends on (backwards). Notice this is the opposite of the regular
     * tallyLooCurrent routine, because now we are trying to tally the
     * incoming instead of the outgoing angular fluxes. */
    if (direction == 1)
    {
        surfID = segment->_mesh_surface_bwd;
        pe_initial = 0;
    }
    else 
    {
        surfID = segment->_mesh_surface_fwd;
        pe_initial = GRP_TIMES_ANG;
    }

    if (surfID != -1)
    {
        /* notice this is polar flux weights, more than just polar weights */
        weights = track->getPolarWeights();
        polar_fluxes = track->getPolarFluxes();

        /* Defines index: instead of 0, 1 now we have 2 (theta less
         * than pi/2), 3 (theta larger than pi/2) */
        if (track->getPhi() < PI / 2.0)
            index = 2;
        else
            index = 3;
        
        /* Obtains the surface that the segment crosses */
        meshSurface = meshSurfaces[surfID];
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
                currents[e] += polar_fluxes[pe] * weights[p] / 2.0;
                pe++;
            }
        }
        meshSurface->incrementCurrents(currents);
    }
    return;
}


/* Effect: accumulates weights for each surface's each
 * quadrature. 0,1,2 are the surface length accumulator (it needs a
 * 1/sin or 1/cos to convert delta_m which is normal distance b/w
 * tracks to the surface lenth). 3,4,5 are weight accumulator (exactly
 * the weight used to accumulate currents).
 */
void Solver::tallyLooWeight(Track *track, segment *segment, 
                             MeshSurface **meshSurfaces, int direction)
{
    int index = 0, opposite_index = 1, surfID;

    /* Get the ID of the surface that the segment ends on (forward), starts
     * on (backwards)*/
    if (direction == 1)
        surfID = segment->_mesh_surface_fwd;
    else 
        surfID = segment->_mesh_surface_bwd;

    if (surfID != -1)
    {
        /* Defines index */
        if (track->getPhi() > PI / 2.0)
            index = 1;    
        opposite_index = 1 - index;
 
        tallyLooWeightSingle(track, segment, meshSurfaces, surfID, index);
        tallyAsTrackedLengthSingle(track, segment, meshSurfaces, surfID, index,
                                   opposite_index);

    }
}

void Solver::tallyLooWeightSingle(Track *track, segment *segment, 
                                  MeshSurface **meshSurfaces,
                                  int surfID, int index)
{        
    int p;
    double wt;
    MeshSurface *meshSurface;

    /* Obtains the surface that the segment crosses */
    meshSurface = meshSurfaces[surfID];

    /* notice this is polar flux weights, more than just polar weights */
    double *weights = track->getPolarWeights();
    double *sinThetaP = _quad->getSinThetas();

    for (p = 0; p < NUM_POLAR_ANGLES; p++)
    {
        wt = 0.5 * weights[p] / TWO_PI / sinThetaP[p];
        meshSurface->incrementTotalWt(wt, index);
    }
}

void Solver::tallyAsTrackedLengthSingle(Track *track, segment *segment, 
                                        MeshSurface **meshSurfaces, 
                                        int surfID, int index, 
                                        int opposite_index)
{        
    int p;
    double cosTheta, sinTheta, wt2;
    MeshSurface *meshSurface;

    /* Obtains the surface that the segment crosses */
    meshSurface = meshSurfaces[surfID];

    /* notice this is polar flux weights, more than just polar weights */
    double *weights = track->getPolarWeights();
    double *sinThetaP = _quad->getSinThetas();

    /* Cell ID */
    int i = surfID / 8;
    int y = i / _cw;
    int x = i % _cw;
    int s = surfID % 8;

    cosTheta = fabs(cos(track->getPhi()));
    sinTheta = fabs(sin(track->getPhi()));

    for (p = 0; p < NUM_POLAR_ANGLES; p++)
    {
        wt2 = 0.5 * weights[p] / TWO_PI / sinThetaP[p];

        if ((s  == 0) || (s == 2))
            meshSurface->incrementAsTrackedLength(wt2 / cosTheta);
        else if (s < 4)
            meshSurface->incrementAsTrackedLength(wt2 / sinTheta);
        else
        {
            _num_crn += 1;
            wt2 *= 0.5;

            if (s < 5)
            {
                meshSurfaces[surfID - 4]
                    ->incrementAsTrackedLength(wt2 / cosTheta);
                meshSurfaces[surfID - 3]
                    ->incrementAsTrackedLength(wt2 / sinTheta);

                if (x > 0)
                {
                    meshSurfaces[(i - 1) * 8 + 1]
                        ->incrementAsTrackedLength(wt2 / sinTheta);
                }
                // FIXME 
                else //if (_geom->getMesh()->getBoundary(0) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 1]->incrementAsTrackedLength
                        (wt2 / sinTheta);
                }

                if (y < _ch -1) 
                {
                    meshSurfaces[(i + _cw) * 8 + 0]
                        ->incrementAsTrackedLength(wt2 / cosTheta);
                }
                else //if (_geom->getMesh()->getBoundary(1) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 0]->incrementAsTrackedLength
                        (wt2 / cosTheta);
                }

            }
            else if (s < 6)
            {
                meshSurfaces[surfID - 4]
                    ->incrementAsTrackedLength(wt2 / sinTheta);
                meshSurfaces[surfID - 3]
                    ->incrementAsTrackedLength(wt2 / cosTheta);

                if (x < _cw - 1)
                {
                    meshSurfaces[(i + 1) * 8 + 1]
                        ->incrementAsTrackedLength(wt2 / sinTheta);
                } 
                else //if (_geom->getMesh()->getBoundary(2) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 1]->incrementAsTrackedLength
                        (wt2 / sinTheta);
                } 
                    
                if (y < _ch -1) 
                {
                    meshSurfaces[(i + _cw) * 8 + 2]
                        ->incrementAsTrackedLength(wt2 / cosTheta);
                }
                else //if (_geom->getMesh()->getBoundary(1) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 2]->incrementAsTrackedLength
                        (wt2 / cosTheta);
                }                      
            }
            else if (s < 7)
            {
                meshSurfaces[surfID - 4]
                    ->incrementAsTrackedLength(wt2 / cosTheta);
                meshSurfaces[surfID - 3]
                    ->incrementAsTrackedLength(wt2 / sinTheta);

                if (x < _cw - 1)
                {
                    meshSurfaces[(i + 1) * 8 + 3]
                        ->incrementAsTrackedLength(wt2 / sinTheta);
                } 
                else //if (_geom->getMesh()->getBoundary(2) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 3]->incrementAsTrackedLength
                        (wt2 / sinTheta);
                }                      

                if (y > 0) 
                {
                    meshSurfaces[(i - _cw) * 8 + 2]
                        ->incrementAsTrackedLength(wt2 / cosTheta);
                }
                else //if (_geom->getMesh()->getBoundary(3) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 2]->incrementAsTrackedLength
                        (wt2 / cosTheta);
                }
            }
            else
            {
                meshSurfaces[surfID - 4]
                    ->incrementAsTrackedLength(wt2 / sinTheta);
                meshSurfaces[surfID - 7]
                    ->incrementAsTrackedLength(wt2 / cosTheta);

                if ((x > 0))// && (y > 0))
                {
                    meshSurfaces[(i - 1) * 8 + 3]
                        ->incrementAsTrackedLength(wt2 / sinTheta);
                }
                else //if (_geom->getMesh()->getBoundary(0) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 3]->incrementAsTrackedLength
                        (wt2 / sinTheta);
                }

                if (y > 0) 
                {
                    meshSurfaces[(i - _cw) * 8 + 0]
                        ->incrementAsTrackedLength(wt2 / cosTheta);
                }
                else //if (_geom->getMesh()->getBoundary(3) == REFLECTIVE)
                {
                    meshSurfaces[i * 8 + 0]->incrementAsTrackedLength
                        (wt2 / cosTheta);
                }                  
            }
        }
    }
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
    FlatSourceRegion* fsr;
    double w, h;
    for (int i = 0; i < _cw * _ch; i++)
    {
        meshCell = _geom->getMesh()->getCells(i);
        for (int s = 8 * i; s < 8 * i + 4; s++)
        {
            meshSurfaces[s]->setTotalWt((meshSurfaces[s]->getTotalWt(0) + 
                                         meshSurfaces[s]->getTotalWt(1)), 2);
        }

        w = meshSurfaces[i * 8 + 1]->getAsTrackedLength();
        h = meshSurfaces[i * 8 + 2]->getAsTrackedLength();
        meshCell->setATWidth(w);
        meshCell->setATHeight(h);
        meshCell->setATL(0.5 * sqrt(w * w + h * h) / P0);

        // Compute mesh cell as-tracked volume by adding up all the FSRs'. 
        double vol = 0;
        for (auto it = meshCell->getFSRs()->begin(); 
             it != meshCell->getFSRs()->end(); it++)
        {
            fsr = &_flat_source_regions[*it];
            vol += fsr->getVolume();
        }
        meshCell->setATVolume(vol);

        if (i == 0)
        {
            log_printf(NORMAL, "physical: %f %f %f; as-tracked: %f %f %f", 
                       meshCell->getWidth(),  meshCell->getHeight(), 
                       meshCell->getVolume(), meshCell->getATWidth(), 
                       meshCell->getATHeight(), meshCell->getATVolume());
        }
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
    int j, k, s, p, e, pe;
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
    if (_run_loo && (!(_first_diffusion && (moc_iter < _num_first_diffusion))))
    {
        tally = 2;
        log_printf(INFO, "tally quadrature currents");
    }
    else if (_run_cmfd || _run_loo)
    {
        tally = 1;
        log_printf(INFO, "tally partial/net currents");
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
        /* Loop over each thread */
        log_printf(ACTIVE, "At the beginning of a sweep");
        for (j = 0; j < _num_azim; j++)
        {
            /* Loop over all tracks for this azimuthal angles */
            for (k = 0; k < _num_tracks[j]; k++) 
            {
                /* Initialize local pointers to important data structures */
                track = &_tracks[j][k];

                track->printOutInfo();

                segments = track->getSegments();
                num_segments = track->getNumSegments();
                weights = track->getPolarWeights();
                polar_fluxes = track->getPolarFluxes();

                /* Store all the incoming angular fluxes */
                //if ((tally == 2) && (!_reflect_outgoing))
                if (tally == 2)
                {
                    tallyLooCurrentIncoming(track, segments.at(0), 
                                            meshSurfaces, 1);
                    tallyLooCurrentIncoming(track, 
                                            segments.at(num_segments-1), 
                                            meshSurfaces, -1);
                }
                /* FIXME: add in tallyCmfdCurrentIncoming */



                /* Loop over each segment in forward direction */
                for (s = 0; s < num_segments; s++) 
                {
                    segment = segments.at(s);
                    fsr = &_flat_source_regions[segment->_region_id];
                    ratios = fsr->getRatios();

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
                                      _pre_factor_array[index + 2 * p +1]));
                            fsr_flux[e] += delta * weights[p];
                            polar_fluxes[pe] -= delta;
                            pe++;
                        }
                    }

#else
                    /* Loop over all polar angles and energy groups */
                    for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
                    {
                        fsr_flux[e] = 0.0;
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

                log_printf(DEBUG, "flux (end of forward) %f", polar_fluxes[0]);

                /* Transfer fluxes to outgoing track, or store them */
                if (_reflect_outgoing)                    
                {
                    track->getTrackOut()->setPolarFluxes(track->isReflOut(),
                                                         0, polar_fluxes);
                }
                else
                    track->setNewFluxes(0, polar_fluxes);                      

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
                        fsr_flux[e] = 0.0;
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
                    log_printf(DEBUG, "flux (end of bacward) %f", 
                               polar_fluxes[0]);

                    /* if segment crosses a surface in bwd direction, 
                       tally quadrature flux for LOO acceleration */
                    if (tally == 2)
                        tallyLooCurrent(track, segment, meshSurfaces, -1);
                    else if (tally == 1)
                        tallyCmfdCurrent(track, segment, meshSurfaces, -1);
						
                    /* Increments the scalar flux for this FSR */
                    fsr->incrementFlux(fsr_flux);
                } /* end of this segment, move on to the next segment */

                /* Transfers fluxes to incoming track, or store them */
                if (_reflect_outgoing)
                {
                    track->getTrackIn()->setPolarFluxes(track->isReflIn(), 
                                                        GRP_TIMES_ANG, 
                                                        polar_fluxes);
                }
                else
                    track->setNewFluxes(GRP_TIMES_ANG, polar_fluxes);

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

            } /* end of a backward segment */
        } /* end of a track (forward & backward segments) */

        log_printf(ACTIVE, "At the end of a sweep, the outgoing fluxes are:");
        if (!_reflect_outgoing)
        {
            /* Loop over each track, updating all the incoming angular fluxes */
            for (j = 0; j < _num_azim; j++) 
            {
                for (k = 0; k < _num_tracks[j]; k++) 
                {
                    track = &_tracks[j][k];

                    /* Transfers outgoing flux to its reflective
                     * track's incoming for the forward direction */
                    track->getTrackOut()
                        ->setPolarFluxes(track->isReflOut(), 0,
                                         track->getNewFluxes());

                    /* Transfers outgoing flux to its reflective
                     * track's incoming for the backward direction */
                    track->getTrackIn()
                        ->setPolarFluxes(track->isReflIn(), GRP_TIMES_ANG,  
                                         track->getNewFluxes());
                    track->printOutNewFluxes();
                    track->zeroNewFluxes();            
                }
            }
        }

        /* If more than one iteration is requested, we only computes source for
         * the last iteration, all previous iterations are considered to be 
         * converging boundary fluxes */
        if (i == max_iterations - 1)
        {
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
                    double f = FOUR_PI * ratios[e] + 
                        (scalar_flux[e] / (2.0 * sigma_t[e] * volume)); 
                    fsr->setFlux(e, f);
                }
            }
		
            normalizeFlux(moc_iter+0.5);

            /* computes new _k_eff; it is important that we compute new k 
             * before computing new source */
            _k_eff = computeKeff(i);
            updateSource(i);
        }
    } /* exit iteration loops */
		
    return;
} /* end of MOCsweep */


void Solver::zeroVacuumBoundaries() {
    Track* track;

    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++) {
            track = &_tracks[i][j];

            track->getTrackOut()
                ->resetPolarFluxes(track->isReflOut(), 0);
            
            track->getTrackIn()
                ->resetPolarFluxes(track->isReflIn(), GRP_TIMES_ANG);
        }
    }
}



/* Normalizes the scalar flux in each FSR and angular flux track to that total 
 * volume-integrated fission source add up to the total volume. */
void Solver::normalizeFlux(double moc_iter)
{			
#if 1
    int counter = 0; 
    double fission_source = 0, factor;
    double *cell_source;

    cell_source = _cmfd->computeCellSourceFromFSR(moc_iter);
    for (int i = 0; i < _cw * _ch; i++)
    {
        if (cell_source[i] > 1e-10)
        {
            counter += 1;
            fission_source += cell_source[i];
        }
    }

    /* Renormalize scalar fluxes in each region */
    factor = counter / fission_source;
    log_printf(ACTIVE, "iter %.1f normalization factor = %f", moc_iter, factor);

    for (int r = 0; r < _num_FSRs; r++)
        _flat_source_regions[r].normalizeFluxes(factor);

    /* Renormalize angular boundary fluxes for each track */
    for (int i = 0; i < _num_azim; i++) 
    {
        for (int j = 0; j < _num_tracks[i]; j++)
            _tracks[i][j].normalizeFluxes(factor);
    }

    /* Renormalize leakage term for each of the four exterior surfaces */
    for (int s = 1; s < 5; s++)
    {
        for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        {
            _geom->getSurface(s)->updateLeakage(factor, e);
        }
    }


    /* This block normalizes tallied current on each surface. To get
     * geometry_2x2.xml with 2 vacuum BC to converge, the tallied
     * current cannot be normalized */
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
                    /* FIXME: I changed 2 to 4, which cases a
                     * SIGFPE, Arithmetic exception */
                    for (int ind = 0; ind < 2; ind++)
                    {
                        meshCell->getMeshSurfaces(s)->updateQuadCurrent(
                            factor, e, ind);
                    }
                }
            }
        }		
    }
#endif
    return;
}

/* Compute the source for each region */
void Solver::updateSource(int moc_iter)
{
#if 1
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

        scalar_flux = fsr->getFlux();
        source = fsr->getSource();
        material = fsr->getMaterial();
        nu_sigma_f = material->getNuSigmaF();
        chi = material->getChi();
        sigma_s = material->getSigmaS();

        /* sigma_s[G * NUM_ENERGY_GROUPS + g] is g to G */
        /*
        log_printf(NORMAL, "sigma_s: 1->1 %f, 1->2 %f, 2->1 %f, 2->2 %f",
                   sigma_s[0], sigma_s[NUM_ENERGY_GROUPS], 
                   sigma_s[1], sigma_s[NUM_ENERGY_GROUPS + 1]);
        */

        /* Compute total fission source for current region */
        fission_source = 0;
        start_index = material->getNuSigmaFStart();
        end_index = material->getNuSigmaFEnd();
        for (int e = start_index; e < end_index; e++)
            fission_source += scalar_flux[e] * nu_sigma_f[e];

        /* Compute total scattering source for group G */
        for (int G = 0; G < NUM_ENERGY_GROUPS; G++) 
        {
            scatter_source = 0;

            start_index = material->getSigmaSStart(G);
            end_index = material->getSigmaSEnd(G);

            for (int g = start_index; g < end_index; g++)
            {
                scatter_source += sigma_s[G * NUM_ENERGY_GROUPS + g]
                    * scalar_flux[g];
            }

            /* Set the total source for region r in group G */
            source[G] = (chi[G] * fission_source / _k_eff + scatter_source) 
                * ONE_OVER_FOUR_PI;

            /* If negative FSR sources show up in the first 5 iterations, set them to zero */
            if ((source[G] < 0) && (moc_iter <= 5))
                source[G] = 0.0;
        }
    }

    /* Update pre-computed source / sigma_t ratios */
    computeRatios();
#endif
}


double Solver::runLoo(int moc_iter)
{
    double loo_keff;

    _k_half = computeKeff(100);

    _cmfd->computeXS();
			 
    if ((moc_iter < _num_first_diffusion) && _first_diffusion)
    {
        _cmfd->computeDs(moc_iter);
        loo_keff = _cmfd->computeCMFDFluxPower(DIFFUSION, moc_iter);
    }
    else 
    {
        _cmfd->computeQuadFlux();

        storeMOCBoundaryFlux();

        _cmfd->computeQuadSrc();
        loo_keff = _cmfd->computeLooFluxPower(moc_iter, _k_eff);
    }

    return loo_keff;
}

double Solver::runCmfd(int moc_iter)
{
    double cmfd_keff;

    _k_half = computeKeff(100);

    /* compute cross sections and diffusion coefficients */
    //_cmfd->computeCurrent();
    _cmfd->computeXS();
    _cmfd->computeDs(moc_iter);

    /* Check for neutron balance */
    if (moc_iter == 10000)
    {
        checkNeutronBalance();
        checkNeutronBalanceWithDs();
    }

    /* Run diffusion problem on initial geometry */
    if ((moc_iter < _num_first_diffusion) && _first_diffusion)
        cmfd_keff = _cmfd->computeCMFDFluxPower(DIFFUSION, moc_iter);
    else
        cmfd_keff = _cmfd->computeCMFDFluxPower(CMFD, moc_iter);

    return cmfd_keff;
}

double Solver::kernel(int max_iterations) {
    log_printf(INFO, "Starting kernel ...");

    double eps_2;

    /* Initial guess */
    _old_k_effs.push(_k_eff);
    _old_eps_2.push(1.0);

    log_printf(INFO, "Starting guess of k_eff = %f", _k_eff);

    /* Check that each FSR has at least one segment crossing it */
    checkTrackSpacing();

    /* Check and print out boundary conditions */
    if ((_run_cmfd || _run_loo) && !(_acc_after_MOC_converge))
        checkBoundary();

    /* Set scalar flux to unity for each region */
    oneFSRFlux();
    initializeTrackFluxes(1.0);
    normalizeFlux(0);
    updateSource(0);
    initializeTrackFluxes(0);//ONE_OVER_FOUR_PI);

    _cmfd->setCellSource(_cmfd->computeCellSourceFromFSR(0));
    _cmfd->printCellSource(0);

    /* Source iteration loop */
    for (int moc_iter = 0; moc_iter < max_iterations; moc_iter++) 
    {
        log_printf(INFO, "Iteration %d: k_eff = %f", moc_iter, _k_eff);

        if (_run_loo)
            _cmfd->storePreMOCMeshSource(_flat_source_regions);
 
        /* Perform one sweep for no acceleration */
        MOCsweep(_boundary_iteration + 1, moc_iter);

        /* Perform acceleration */
        if ((_run_cmfd) && !(_acc_after_MOC_converge))
            _acc_k = runCmfd(moc_iter);
        else if ((_run_loo) && !(_acc_after_MOC_converge))
            _acc_k = runLoo(moc_iter);

        _cmfd->incrementMOCIter();

        /* Update FSR's scalar flux and boundary angular fluxes */
        if ((_run_cmfd || _run_loo) && !(_acc_after_MOC_converge))
        {
            prolongation(moc_iter);
            normalizeFlux(moc_iter+1);
            _cmfd->printCellSource(moc_iter+1);
            _k_eff = _acc_k;
            updateSource(moc_iter);
        }

        /* We only store $k^{(m+1)}$; other intermediate keff does not matter */
        _old_k_effs.push(_k_eff);
        if (_old_k_effs.size() == NUM_KEFFS_TRACKED)
            _old_k_effs.pop();

        eps_2 = _cmfd->computeCellSourceNorm();

        _old_eps_2.push(eps_2);
        if (_old_eps_2.size() == NUM_KEFFS_TRACKED)
            _old_eps_2.pop();

        /* Prints out keff & eps, may update keff too based on _update_keff */
        /* FIXME: the following line generates a bad read in
         * Valgrine */
        printToScreen(moc_iter + 1);
        printToLog(moc_iter + 1);
        if (_plot_flux)
            plotFluxes(moc_iter + 1);

        if (moc_iter < 35)
            _cmfd->printCellSource(moc_iter + 1);


        /* Alternative: if (_cmfd->getL2Norm() < _moc_conv_thresh) */
        if (eps_2 < _moc_conv_thresh) 
        {

            _cmfd->printCellSource(1000);

            
            /* Run one steps of acceleration if it is requested to do so */
            if (_acc_after_MOC_converge)
            {
                if (_run_loo)
                    _acc_k = runLoo(10000);
                if (_run_cmfd)
                    _acc_k = runCmfd(10000);
            }

            printToMinimumLog(moc_iter + 1);
            log_printf(NORMAL, "Printed log file to %s", 
                       _log_file.c_str());

            plotEverything(moc_iter);

            for (int j = 0; j < _num_azim; j++)
            {
                for (int k = 0; k < _num_tracks[j]; k++)
                {
                    log_printf(DEBUG, "final boundary flux = %f %f",
                               _tracks[0][0].getPolarFluxes()[0], 
                               _tracks[0][0].getPolarFluxes()[1]);
                }
            }

            return _k_eff;
        }
        /* FIXME: this is an arbitrary number */
        else if (eps_2 > 1e3) 
        {
            printToMinimumLog(moc_iter);
            plotEverything(moc_iter);

            log_printf(NORMAL, "Exit OpenMOC: L2 norm blows up."
                       "  Something bad happens. May need damping.");
            return _k_eff;
        }
    } /* end of moc iterations */

    log_printf(WARNING, "Unable to converge the source after %d iterations",
               max_iterations);

    return _k_eff;
}

void Solver::plotEverything(int moc_iter)
{
    /* plot pin powers */
    if (_compute_powers)
        plotPinPowers();

    _plotter->plotFSRs(_geom->getMesh(), _num_FSRs);

    /* plot CMFD flux and xs */
    if (_run_cmfd && _plotter->plotCurrent() )
    {
        _cmfd->computeXS();
        _cmfd->computeDs(moc_iter);
        _plotter->plotDHats(_geom->getMesh(), moc_iter);
        _plotter->plotXS(_geom->getMesh(), moc_iter);
    }

    /* plot LOO flux and xs */
    if ((_run_loo) && (_plot_loo))
    {
        _plotter->plotQuadFlux(_geom->getMesh(), moc_iter);
        _plotter->plotNetCurrents(_geom->getMesh(), moc_iter);
        _plotter->plotXS(_geom->getMesh(), moc_iter);
    }

    /* plot FSR scalar flux */
    if (_plot_flux)
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


void Solver::printToLog(int moc_iter)
{
    std::ofstream logfile;

    if (moc_iter == 1)
    {
        logfile.open(_log_file.c_str(), std::fstream::trunc);
        logfile << "# iteration,"
                << " cell l2 norm (m, m+1),"
                << " cell l2 norm ratio"
                << " keff relative change,"
                << " #lo iterations, "
                << " keff"
                << std::endl;
    }
    else
    {
        logfile.open(_log_file.c_str(), std::ios::app);
        logfile << moc_iter 
                << " " << _old_eps_2.back()
                << " " << _old_eps_2.front() / _old_eps_2.back()
                << " " << 1.0 -  _old_k_effs.front() / _old_k_effs.back()
                << " " << _cmfd->getNumIterToConv() 
                << " " << std::setprecision(11) << _old_k_effs.back()
                << std::endl;
    }

    logfile.close();
        
}



/*
 * check neutron balance in each mesh cell assuming reflective outter
 * boundary condition FIXME: need to update to include vacuum BC. 
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
            /* get mesh cell */
            meshCell = mesh->getCells(y * cell_width + x);

            for (int e = 0; e < ng; e++)
            {
                leak = 0;
                absorb = 0;
                fis = 0;
                src = 0;

                flux = meshCell->getOldFlux()[e];
                vol = meshCell->getVolume();
				
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
                log_printf(DEBUG, "CMFD cell %d energy %d, residual %.10f"
                           " leak: %e, absorb: %.10f, src: %.10f,"
                           " fis: %.10f, keff: %.10f", 
                           y * cell_width + x, e, residual,
                           leak, absorb, src, fis, fis / (leak + absorb));
                tot_leak += leak;
                tot_absorb += absorb;
                tot_fis += fis;
                tot_src += src;
            }	
        }
    }

    log_printf(INFO, "Assume no boundary condition, keff: %.10f"
               " fis: %f, absorb: %f, leak: %e, src = %f, res = %f", 
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
        log_printf(DEBUG, "CMFD residual %.10f"
                   " fis: %.10f, absorb: %.10f, leak: %.10f, keff: %.10f", 
                   residual, fis, absorb, leak, fis / (leak + absorb));
    } /* end of looping through cells */

    return;
}


FlatSourceRegion* Solver::getFSRs(){
    return _flat_source_regions;
}

int Solver::getNumFSRs(){
    return _num_FSRs;
}
