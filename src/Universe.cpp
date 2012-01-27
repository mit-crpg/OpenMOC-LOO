/*
 * Universe.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Universe.h"


/* _n keeps track of the number universes of instantiated */
int Universe::_n = 0;

/**
 * Universe constructor
 * @param id the universe id
 */
Universe::Universe(const int id) {
	_uid = _n;
	_id = id;
	_n++;
}


/**
 * Destructor
 */
Universe::~Universe() {
	_cells.clear();
}


/**
 * Adds a cell to this universe
 * @param cell the cell id
 */
void Universe::addCell(Cell* cell) {
	try {
		_cells.push_back(cell);
		log_printf(INFO, "Added cell with id = %d to universe with id = %d\n",
				cell->getId(), _id);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add cell with id = %d to universe with"
				" id = %d. Backtrace:\n%s", cell, _id, e.what());
	}
}


/**
 * Return the vector of cells in this universe
 * @return vector of cell ids
 */
std::vector<Cell*> Universe::getCells() const {
    return _cells;
}


/**
 * Returns the universe's uid
 * @return the universe's uid
 */
int Universe::getUid() const {
	return _uid;
}


/**
 * Return the id for this universe
 * @return the universe id
 */
int Universe::getId() const {
    return _id;
}


/**
 * Return the number of cells in this universe
 * @return the number of cells
 */
int Universe::getNumCells() const {
    return _cells.size();
}


/**
 * Return a pointer to the origin for this cell (in global coordinates)
 * @return the origin of the cell
 */
Point* Universe::getOrigin() {
    return &_origin;
}


/**
 * Set the origin for this universe
 * @param origin the origin point
 */
void Universe::setOrigin(Point* origin) {
    _origin.setX(origin->getX());
    _origin.setY(origin->getY());
}



/**
 * Convert the member attributes of this universe to a character array
 * @param a character array reprsenting the universe
 */
std::string Universe::toString() {
	std::stringstream string;

	string << "Universe id = " << _id << ", num cells = " << _cells.size() <<
			", cell ids = ";

	for (int c = 0; c < (int)_cells.size(); c++)
		string << _cells.at(c)->getId() << ", ";

	string << "\n";

	return string.str();
}
