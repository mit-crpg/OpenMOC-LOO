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
Universe::Universe(int id) {
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
void Universe::addCell(int cell) {
	_cells.push_back(cell);
	_num_cells++;
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
 * Return the level for this universe (base level = 0)
 * @return the universe level
 */
int Universe::getLevel() const {
    return _level;
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
 * Set the level of this universe
 * @param level the universe level (base = 0)
 */
void Universe::setLevel(int level) {
    _level = level;
}


/**
 * Set the origin for this universe
 * @param origin the origin point
 */
void Universe::setOrigin(Point* origin)
{
    _origin.setX(origin->getX());
    _origin.setY(origin->getY());
}
