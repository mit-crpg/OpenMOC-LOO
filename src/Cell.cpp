/*
 * Cell.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Cell.h"

/* _n keeps track of the number cells of instantiated */
int Cell::_n = 0;

/**
 * Default Cell constructor
 * @param id the cell id
 * @param type the type of cell
 * @param universe the universe this cell is in
 * @param parent_cell the cell within which this cell resides
 * @param num_surfaces the number of surfaces in this cell
 * @param surfaces the surface id
 */
Cell::Cell(int id, cellType type, int universe, int num_surfaces, std::vector<int> surfaces) {
	_uid = _n;
	_id = id;
	_type = type;
	_universe = universe;
	_num_surfaces = num_surfaces;
	_surfaces = surfaces;
	_n++;
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
		log_printf(ERROR, 
			   "Unable to add surface with id = %e to cell"
			   "with id = %d. Backtrace:\n%s\n", 
			   surface, _id, e.what());
	}
}


/**
 * Return the cell's uid
 * @return the cell's uid
 */
int Cell::getUid() const {
	return _uid;
}



/**
 * Return the cell's id
 * @return the cell's id
 */
int Cell::getId() const {
	return _id;
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
 * Return the number of surfaces in the cell
 * @return the number of surfaces
 */
int Cell::getNumSurfaces() const {
	return _num_surfaces;
}

/**
 * Return the vector of surfaces in the cell
 * @return vector of surface ids
 */
std::vector<int> Cell::getSurfaces() const {
	return _surfaces;
}


/**
 * Set the universe that this cell is inside of
 * @param the universe's id
 */
void Cell::setUniverse(int universe) {
	_universe = universe;
}

/**
 *  CellBasic constructor
 *  @param material the material used to fill this cell
 */
CellBasic::CellBasic(int id, int universe, int num_surfaces, 
		   std::vector<int> surfaces, 
		   int material): 
	Cell(id, MATERIAL, universe, num_surfaces, surfaces) {
	_material = material;
}

/**
 * Return the material in the cell
 * @return the material's id
 */
int CellBasic::getMaterial() const {
	return _material;
}


/**
 * Adjusts the cell's id of the universe this cell is in, the id of the
 * universe that is filling this cell, and the id of the surfaces inside
 * this cell to be the uids of each rather than the ids defined by the
 * input file
 */
void CellBasic::adjustKeys(int universe, int material,
		std::map<int, Surface*> surfaces) {
	_universe = universe;
	_material = material;

	std::vector<int> adjusted_surfaces;

	try {
		/* Adjust surface ids to be the surface uids */
		for (int s = 0; s < (int)_surfaces.size(); s++)
			adjusted_surfaces.push_back(surfaces.at(s)->getUid());

		_surfaces.clear();
		_surfaces = adjusted_surfaces;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust the surface ids to uids for cell "
				"id = %d. Backtrace:\n%s", _id, e.what());
	}
}


/**
 * Convert this cell's attributes to a string format
 * @return a character array of this cell's attributes
 */
const char* CellBasic::toString() {
	std::stringstream string;
	string << "Cell id = " << _id << ", type = MATERIAL, material id = " <<
			_material << ", universe = " << _universe << ", num_surfaces = "
			<< _num_surfaces << " surface ids = ";

	for (int s = 0; s < _num_surfaces; s++)
		string << _surfaces.at(s) << ", ";

	string << std::endl;

	return string.str().c_str();
}


/**
 *  CellFill constructor
 *  @param universe_fill the universe used to fill this cell
 */
CellFill::CellFill(int id, int universe, int num_surfaces, 
		   std::vector<int> surfaces, 
		   int universe_fill): 
	Cell(id, FILL, universe, num_surfaces, surfaces) {
	_universe_fill = universe_fill;
}

/**
 * Return the universe filling this cell
 * @return the universe's id
 */
int CellFill::getUniverseFill() const {
	return _universe_fill;
}

/**
 * Set the universe filling this cell
 * @param the universe's id
 */
void CellFill::setUniverseFill(int universe_fill) {
	_universe_fill = universe_fill;
}


/**
 * Adjusts the cell's id of the universe this cell is in, the id of the
 * universe that is filling this cell, and the id of the surfaces inside
 * this cell to be the uids of each rather than the ids defined by the
 * input file
 */
void CellFill::adjustKeys(int universe, int universe_fill,
		std::map<int, Surface*> surfaces) {
	_universe = universe;
	_universe_fill = universe_fill;

	std::vector<int> adjusted_surfaces;

	try {
		/* Adjust surface ids to be the surface uids */
		for (int s = 0; s < (int)_surfaces.size(); s++)
			adjusted_surfaces.push_back(surfaces.at(s)->getUid());

		_surfaces.clear();
		_surfaces = adjusted_surfaces;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust the surface ids to uids for cell "
				"id = %d. Backtrace:\n%s", _id, e.what());
	}

}


/**
 * Convert this cell's attributes to a string format
 * @return a character array of this cell's attributes
 */
const char* CellFill::toString() {
	std::stringstream string;
	string << "Cell id = " << _id << ", type = FILL, universe_fill = " <<
			_universe_fill << ", universe = " << _universe << ", num_surfaces = "
			<< _num_surfaces;

	string << ", surface ids: ";
	for (int s = 0; s < _num_surfaces; s++)
		string << _surfaces.at(s) << ", ";

	string << "\n";

	return string.str().c_str();
}


