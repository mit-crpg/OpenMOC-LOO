/*
 * LocalCoords.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "LocalCoords.h"

/**
 * LocalCoords constructor
 * @param x the x-coordinate
 * @param y the y-coordinate
 */
LocalCoords::LocalCoords(double x, double y) {
	_coords.setCoords(x, y);
	_next = NULL;
	_prev = NULL;
}


/**
 * LocalCoords destructor recursively deletes all LocalCoords objects
 * behind this one in the linked list
 */
LocalCoords::~LocalCoords() {
	delete _next;
}


/**
 * Return the level (CELL or LATTICE) of this localcoords
 * @return the level (CELL or LATTICE)
 */
coordType LocalCoords::getType() {
	return _type;
}


/**
 * Return the universe id that this coordinate is in
 * @return the universe id
 */
int LocalCoords::getUniverse() const {
    return _universe;
}


/**
 * Return the cell id for this localcoords
 * @return the cell id
 */
int LocalCoords::getCell() const {
	return _cell;
}


/**
 * Return the lattice id that this coordinate is in
 * @return the lattice id
 */
int LocalCoords::getLattice() const {
    return _lattice;
}

/**
 * Return lattice cell along the x-axis that this coordinate is in
 * @return lattice cell x
 */
int LocalCoords::getLatticeX() const {
    return _lattice_x;
}


/**
 * Return lattice cell along the y-axis that this coordinate is in
 * @return lattice cell y
 */
int LocalCoords::getLatticeY() const {
    return _lattice_y;
}


/**
 * Returns the x-coordinate
 * @return the x-coordinate
 */
double LocalCoords::getX() const {
    return _coords.getX();
}


/**
 * Returns the y-coordinate
 * @return the y-coordinate
 */
double LocalCoords::getY() const {
    return _coords.getY();
}


/**
 * Returns a pointer to the coordinates for this localcoord
 * @return pointer the coordinates
 */
Point* LocalCoords::getPoint() {
	return &_coords;
}


/**
 * Return a pointer to the localcoord at the next level if one exists
 * @return pointer to the next localcoord
 */
LocalCoords* LocalCoords::getNext() const {
    return _next;
}


LocalCoords* LocalCoords::getPrev() const {
	return _prev;
}

/**
 * Set the level for this localcoords
 * @param level the level for this localcoords (LATTICE or CELL)
 */
void LocalCoords::setType(coordType type) {
	_type = type;
}


/**
 * Sets the universe id that this coordinate is in
 * @param universe the universe id
 */
void LocalCoords::setUniverse(int universe) {
    _universe = universe;
}


/**
 * Set the cell id for this localcoords
 * @param cell the cell id
 */
void LocalCoords::setCell(int cell) {
	_cell = cell;
}


/**
 * Sets the lattice id that this coordinate is in
 */
void LocalCoords::setLattice(int lattice) {
    _lattice = lattice;
}


/**
 * Sets the lattice cell along the x-axis that this coordinate is in
 * @param lattice_x the x lattice cell
 */
void LocalCoords::setLatticeX(int lattice_x) {
    _lattice_x = lattice_x;
}


/**
 * Sets the lattice cell along the y-axis that this coordinate is in
 * @param lattice_y the y lattice cell
 */
void LocalCoords::setLatticeY(int lattice_y) {
    _lattice_y = lattice_y;
}


/**
 * Set the x-coordinate
 * @param x the x-coordinate
 */
void LocalCoords::setX(double x) {
	_coords.setX(x);
}


/**
 * Set the y-coordinate
 * @param y the y-coordinate
 */
void LocalCoords::setY(double y) {
	_coords.setY(y);
}


/**
 * Sets the pointer to the next localcoord object
 * @param next
 */
void LocalCoords::setNext(LocalCoords* next) {
    _next = next;
}

void LocalCoords::setPrev(LocalCoords* prev) {
	_prev = prev;
}


void LocalCoords::adjustCoords(double delta_x, double delta_y) {

	/* Forward direction along linked list */
	LocalCoords* curr = this;
	while (curr != NULL) {
		curr->setX(curr->getX() + delta_x);
		curr->setY(curr->getY() + delta_y);
		curr = curr->getNext();
	}

	/* Reverse direction along linked list */
	curr = _prev;
	while (curr != NULL) {
		curr->setX(curr->getX() + delta_x);
		curr->setY(curr->getY() + delta_y);
		curr = curr->getPrev();
	}
}


/**
 * Converts this localcoords's attributes to a character array
 * representation
 * @param a character array of its member's attributes
 */
std::string LocalCoords::toString() {

	std::stringstream string;
	LocalCoords* curr = this;

	while (curr != NULL) {
		string << "LocalCoords: level = ";

		if (curr->getType() == UNIV) {
			string << " UNIVERSE, x = " << curr->getX() << ", y = " << curr->getY()
					<< ", universe = " << curr->getUniverse() << ", cell = " <<
					curr->getCell();
		}
		else if (curr->getType() == LAT){
			string << " LATTICE, x = " << curr->getX() << ", y = " << curr->getY()
					<< ", universe = " << curr->getUniverse() << ", lattice = " <<
					curr->getLattice() << ", lattice_x = " << curr->getLatticeX()
					<< ", lattice_y = " << curr->getLatticeY();
		}
		else {
			string << " NONE, x = " << curr->getX() << ", y = " << curr->getY()
					<< ", universe = " << curr->getUniverse() << ", lattice = " <<
					curr->getLattice() << ", lattice_x = " << curr->getLatticeX()
					<< ", lattice_y = " << curr->getLatticeY()
					<< ", cell = " << curr->getCell();
		}

		string << ", next:\n";
		curr = curr->getNext();
	}


	return string.str();
}
