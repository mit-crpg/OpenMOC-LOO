/*
 * Universe.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Universe.h"


/**
 * Universe constructor
 * @param id the universe id
 */
Universe::Universe(const int id) {
	_id = id;
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
void Universe::addCell(const int cell) {
	try {
		_cells.push_back(cell);
		_num_cells++;
		log_printf(INFO, "Added cell with id = %d to universe with id = %d", cell, _id);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add cell with id = %d to universe with id = %d. "
				"Backtrace:\n%s", cell, _id, e.what());
	}
}


/**
 * Return the vector of cells in this universe
 * @return vector of cell ids
 */
std::vector<int> Universe::getCells() const {
    return _cells;
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
    return _num_cells;
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
const char* Universe::toString() {
	std::stringstream string;

	string << "Universe id = " << _id << ", num cells = " << _num_cells <<
			", cell ids = ";

	for (int c = 0; c < _num_cells; c++)
		string << _cells.at(c) << ", ";

	string << std::endl;

	return string.str().c_str();
}
