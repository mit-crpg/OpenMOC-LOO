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
Cell::Cell(int id, cellType type, int universe, int num_surfaces,
		int *surfaces) {

	_uid = _n;
	_id = id;
	_type = type;
	_universe = universe;
	_n++;
	
	/* This empty surface pointer is just a null value for the _surfaces
	 * map. The Geometry will register the actual surface pointer when
	 * the cell is added to the geometry
	 */
	Surface* empty_surface_pointer;
	for (int i = 0; i < num_surfaces; i++)
		_surfaces.insert(std::pair<int, Surface*>(surfaces[i],
							empty_surface_pointer));
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
void Cell::setSurfacePointer(Surface* surface) {
	/* iF _surfaces does not contain this surface id throw an error */
	if (_surfaces.find(surface->getId()) == _surfaces.end() &&
			_surfaces.find(-surface->getId()) == _surfaces.end())

		log_printf(WARNING, "Unable to set surface pointer for cell id = %d "
				"for surface id = %d since cell does not contain this surface\n",
				_id, surface->getId());

	try{
		/* If the cell contains the negative side of the surface */
		if (_surfaces.find(surface->getId()) != _surfaces.end())
			_surfaces[surface->getId()] = surface;

		/* If the cell contains the negative side of the surface */
		else
			_surfaces[-1*surface->getId()] = surface;

		log_printf(INFO, "Set the surface pointer for cell id = %d for "
				"surface id = %d\n", _id, surface->getId());
	}
	catch (std::exception &e) {
		log_printf(ERROR, 
			   "Unable to add surface with id = %d to cell with id = %d. "
			   "Backtrace:\n%s", surface, _id, e.what());
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
	return this->_surfaces.size();
}


/**
 * Return the vector of surfaces in the cell
 * @return vector of surface ids
 */
std::map<int,Surface*> Cell::getSurfaces() const {
	return _surfaces;
}


/**
 * Set the universe that this cell is inside of
 * @param the universe's id
 */
void Cell::setUniverse(int universe) {
	_universe = universe;
}


/*
 * Determines whether a point is contained inside a cell. Queries each surface
 * inside the cell to determine if the particle is on the same side of the
 * surface. This particle is only inside the cell if it is on the same side of
 * every surface in the cell.
 * @param point a pointer to a point
 */
bool Cell::cellContains(Point* point) {

//	std::vector<int> cell_surfaces = cell->getSurfaces();
	std::map<int, Surface*>::iterator iter;

	for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
		if (iter->second->evaluate(point) * iter->first < ON_SURFACE_NEG)
			return false;
	}

	return false;
}


/*
 * Determines whether a point is contained inside a cell. Queries each surface
 * inside the cell to determine if the particle is on the same side of the
 * surface. This particle is only inside the cell if it is on the same side of
 * every surface in the cell.
 * @param point a pointer to a localcoord
 */
bool Cell::cellContains(LocalCoords* coords) {
	return this->cellContains(coords->getPoint());
}


/**
 *  CellBasic constructor
 *  @param material the material used to fill this cell
 */
CellBasic::CellBasic(int id, int universe, int num_surfaces, 
		   int *surfaces, int material):
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
void CellBasic::adjustKeys(int universe, int material) {

	_universe = universe;
	_material = material;
	std::map<int, Surface*> adjusted_surfaces;

	try {
		/* Adjust surface ids to be the positive/negative surface uids */
		std::map<int, Surface*>::iterator iter;
		for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
			/* Positive side of surface, use positive uid */
			if (iter->first > 0)
				adjusted_surfaces.insert(std::pair<int, Surface*>
								(iter->second->getUid(), iter->second));
			/* Negative side of surface, use negative uid */
			else
				adjusted_surfaces.insert(std::pair<int, Surface*>
								(iter->second->getUid()*-1, iter->second));
		}

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
std::string CellBasic::toString() {

	std::stringstream string;

	string << "Cell id = " << _id << ", type = MATERIAL, material id = " <<
			_material << ", universe = " << _universe << ", num_surfaces = "
	       << getNumSurfaces() << ", surface ids = ";

	std::map<int, Surface*>::iterator iter;
	for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
		string << iter->first << ", ";

	string << "\n";

	return string.str();
}


/**
 *  CellFill constructor
 *  @param universe_fill the universe used to fill this cell
 */
CellFill::CellFill(int id, int universe, int num_surfaces, 
		   int *surfaces, int universe_fill):
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
void CellFill::adjustKeys(int universe, int universe_fill) {

	_universe = universe;
	_universe_fill = universe_fill;

	std::map<int, Surface*> adjusted_surfaces;

	try {
		/* Adjust surface ids to be the positive/negative surface uids */
		std::map<int, Surface*>::iterator iter;
		for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
			/* Positive side of surface, use positive uid */
			if (iter->first > 0)
				adjusted_surfaces.insert(std::pair<int, Surface*>
									(iter->second->getUid(), iter->second));
			/* Negative side of surface, use negative uid */
			else
				adjusted_surfaces.insert(std::pair<int, Surface*>
									(iter->second->getUid()*-1, iter->second));
		}

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
std::string CellFill::toString() {

	std::stringstream string;

	string << "Cell id = " << _id << ", type = FILL, universe_fill = " <<
			_universe_fill << ", universe = " << _universe << ", num_surfaces = "
	       << getNumSurfaces();

	std::map<int, Surface*>::iterator iter;
	string << ", surface ids = ";
	for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
		string << iter->first << ", ";

	string << "\n";

	return string.str();
}


