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
FlatSourceRegion::FlatSourceRegion() { }


/**
 * Default destructor
 */
FlatSourceRegion::~FlatSourceRegion() { }

/**
 * Returns this region's id. This id must correspond to the id computed
 * using a region index schema based on the layout of the geometry
 * @return the region's id
 */
int FlatSourceRegion::getId() const {
    return _id;
}


/**
 * Returns the id of the cell that this flat source region is in
 * @return the cell id
 */
int FlatSourceRegion::getCellId() const {
    return _cell_id;
}


/**
 * Returns the id of the material that this flat source region has in it
 * @return the material id
 */
int FlatSourceRegion::getMaterialId() const {
    return _material_id;
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
 * Returns an array of the multi-energy group fluxes tallied in this region
 * during the previous iteration in the solver
 * @return a flux array
 */
double* FlatSourceRegion::getOldFlux() {
    return _old_flux;
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
 * Returns a array of the multi-energy group source term in this region from
 * the previous iteration computed by the solver
 * @return the old source array
 */
double* FlatSourceRegion::getOldSource() {
    return _old_source;
}


/**
 * Sets this region's cell id
 * @param cell_id this region's cell id
 */
void FlatSourceRegion::setCellId(int cell_id) {
    _cell_id = cell_id;
}


/**
 * Sets this region's material id
 * @param material_id the id of the material inside this region
 */
void FlatSourceRegion::setMaterialId(int material_id) {
    _material_id = material_id;
}


/**
 * Sets this region's id
 * @param id the region id
 */
void FlatSourceRegion::setId(int id){
    _id = id;
}


/**
 * Sets the reigon's volume
 * @param volume the region's volume
 */
void FlatSourceRegion::setVolume(double volume) {
    _volume = volume;
}
