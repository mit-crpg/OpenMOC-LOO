/*
 * Cell.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Cell.h"


/**
 * Cell constructor
 * @param id the cell id
 * @param type the type of cell
 * @param universe_fill the id of the universe filling this cell
 * @param material the material filling this cell (0 if filled by universe)
 * @param num_surfaces the number of surfaces in this cell
 */
Cell::Cell(int id, cellType type, int num_surfaces) {
	_id = id;
	_type = type;
	_num_surfaces = num_surfaces;

	/* The following are default values. These values are dependent on
	 * whether the cell is of MATERIAL or FILL type
	 */
	_universe_fill = -1E5;
	_universe = -1E5;
	_material = -10E5;
}


/**
 * Destructor frees all surfaces making up cell
 */
Cell::~Cell() {
	_surfaces.clear();
}


/**
 * Add a surface to the cell
 * @param surface the surface id
 */
void Cell::addSurface(int surface) {
	try{
		_surfaces.push_back(surface);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add surface with id = %e to cell with id = %d. Backtrace:"
				"\n%s\n", surface, _id, e.what());
	}
}


/**
 * Return the cell's id
 * @return the cell's id
 */
int Cell::getId() const {
    return _id;
}


/**
 * Return the material in the cell
 * @return the material's id
 */
int Cell::getMaterial() const {
    return _material;
}


/**
 * Return the number of surfaces in the cell
 * @return the number of surfaces
 */
int Cell::getNumSurfaces() const {
    return _num_surfaces;
}


/**
 * Return the cell' parent
 * @return the parent cell id
 */
int Cell::getParentCell() const {
    return _parent_cell;
}


/**
 * Return the vector of surfaces in the cell
 * @return vector of surface ids
 */
std::vector<int> Cell::getSurfaces() const {
    return _surfaces;
}


/**
 * Return the cell type (FILL or MATERIAL)
 * @return the cell type
 */
cellType Cell::getType() const {
    return _type;
}


/**
 * Return the universe that this cell is in
 * @return the universe id
 */
int Cell::getUniverse() const {
    return _universe;
}

/**
 * Return the universe which is filling this cell
 * @param the universe id
 */
int Cell::getUniverseFill() const {
    return _universe_fill;
}


/**
 * Sets the material inside the cell
 * @param material the material's id
 */
void Cell::setMaterial(int material) {
    _material = material;
}


/**
 * Set the universe that this cell is inside of
 * @param the universe's id
 */
void Cell::setUniverse(int universe) {
    _universe = universe;
}


/**
 * Sets the universe which fills this cell (if cell is of FILL type only)
 * @param the universe's id
 */
void Cell::setUniverseFill(int universeFill) {
    _universe_fill = universeFill;
}


/**
 * Set the parent cell for this cell
 * @param parentCell the parent cell id
 */
void Cell::setParentCell(int parentCell) {
    _parent_cell = parentCell;
}
