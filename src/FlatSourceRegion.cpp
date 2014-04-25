/*
 * FlatSourceRegion.cpp
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "FlatSourceRegion.h"

/**
 * FlatSourceRegion constructor
 */
FlatSourceRegion::FlatSourceRegion() 
{
    _material = NULL;
    _volume = 0.0;
    _quad_id = -1;

    /* Initializes region's other attributes */
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++) 
    {
        _flux[e] = 0.0;
        _source[e] = 0.0;
    }

#if USE_OPENMP
    omp_init_lock(&_flux_lock);
#endif
}


/**
 * Default destructor
 */
FlatSourceRegion::~FlatSourceRegion() {
#if USE_OPENMP
    omp_destroy_lock(&_flux_lock);
#endif
}

/**
 * Returns this region's id. This id must correspond to the id computed
 * using a region index schema based on the layout of the geometry
 * @return the region's id
 */
int FlatSourceRegion::getId() const {
    return _id;
}

int FlatSourceRegion::getQuadId() const {
    return _quad_id;
}

/**
 * Returns this mesh cell's id that this FSR belongs to. 
 * @return the mesh cell's id
 */
int FlatSourceRegion::getMeshCellId() const {
    return _mesh_cell_id;
}

/**
 * Returns a pointer to this region's material
 * @return a pointer to a material
 */
Material* FlatSourceRegion::getMaterial() {
    return _material;
}


/**
 * Returns the region's volume as computed by the number of track segments
 * and their cumulative lengths in this region
 * @preturn the region's volume
 */
double FlatSourceRegion::getVolume() const {
    return _volume;
}


/**
 * Returns an array of the multi energy group fluxes tallied in this region
 * @return a flux array
 */
double* FlatSourceRegion::getFlux() {
    return _flux;
}

/**
 * Returns a array of the multi-energy group source term in this region
 * computed by the solver
 * @return the old source array
 */
double* FlatSourceRegion::getSource() {
    return _source;
}

/**
 * Return an array of ratios of source / sigma_t for this flat source region
 * for each energy group
 * @return array of source / sigma_t ratio
 */
double* FlatSourceRegion::getRatios() {
    return _ratios;
}


/**
 * Sets this region's id
 * @param id the region id
 */
void FlatSourceRegion::setId(int id){
    _id = id;
}

void FlatSourceRegion::setQuadId(int id){
    _quad_id = id;
}

/**
 * Sets the id for the mesh cell that this FSR belongs to 
 * @param id the mesh cell id
 */
void FlatSourceRegion::setMeshCellId(int id){
    _mesh_cell_id = id;
}

/**
 * Set this region's material
 * @param material pointer to a material
 */
void FlatSourceRegion::setMaterial(Material* material) {
    _material = material;
}


/**
 * Sets the reigon's volume
 * @param volume the region's volume
 */
void FlatSourceRegion::setVolume(double volume) {
    _volume = volume;
}


/**
 * Increment this FSR's volume by some amount corresponding to a segment length
 * @param volume the amount to increment by
 */
void FlatSourceRegion::incrementVolume(double volume) {
    assert(volume > 0);
    _volume += volume;
}


/**
 * Set the scalar flux for one energy group inside this FSR
 * @param energy the energy group index
 * @param flux the scalar flux
 */
void FlatSourceRegion::setFlux(int energy, double flux) {
    assert(energy > -1);
    assert(energy < NUM_ENERGY_GROUPS);
    assert(flux > - 1e-10);
    _flux[energy] = flux;
    return;
}


/**
 * Increment the scalar flux for one energy group inside this FSR
 * @param energy the energy group index
 * @param flux the scalar flux
 */
void FlatSourceRegion::incrementFlux(int energy, double flux) {
    assert(energy > -1);
    assert(energy < NUM_ENERGY_GROUPS);
    assert(flux > -1e-10);

    _flux[energy] += flux;
    return;
}


/**
 * Increment the scalar flux for all energy groups inside this FSR
 * @param flux the scalar flux for all energy groups
 */
void FlatSourceRegion::incrementFlux(double* flux) {
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        _flux[e] += flux[e];

    return;
}

/**
 * Set the source for one energy group for this FSR
 * @param energy the energy group index
 * @param source the source
 */
void FlatSourceRegion::setSource(int energy, double source) {
    assert(energy > -1);
    assert(energy < NUM_ENERGY_GROUPS);
    assert(source > -1e-10);
    _source[energy] = source;
    return;
}

/**
 * Normalizes all of the scalar flux values by multiplying by a factor
 * @param factor the factor to scale the flux by
 */
void FlatSourceRegion::normalizeFluxes(double factor) {

    /* Loop over all energy groups */
    for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
        _flux[e] *= factor;

    return;
}


/**
 * Compute source / sigma_t for each energy group
 */
void FlatSourceRegion::computeRatios() {

    double* sigma_t = _material->getSigmaT();

    for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
        _ratios[e] = _source[e] / sigma_t[e];
    }

    return;
}


/**
 * Compute the volumetric fission rate in this flat source region by adding
 * up the fission rates in each energy group. This method assumes that fixed
 * source iteration has already been run since it uses the flux stored in
 * this region
 */
double FlatSourceRegion::computeFissionRate() {
    double power = 0.0;
    double* sigma_f = _material->getNuSigmaF();

    /* Add the fission rates from each energy group */
    for (int e=0; e < NUM_ENERGY_GROUPS; e++)
        power += sigma_f[e] * _flux[e];

    /* Multiply by volume of FSR */
    power *= _volume;

    return power;
}

/*
void FlatSourceRegion::setBoundaryUpdate(int group, int ind, double bu)
{
    _boundary_update[group * 2 + ind] = bu;
}

double FlatSourceRegion::getBoundaryUpdate(int group, int ind)
{
    return _boundary_update[group * 2 + ind];
}
*/
