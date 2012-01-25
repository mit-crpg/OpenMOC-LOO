/*
 * Lattice.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Lattice.h"


/**
 * Lattice constructor
 * @param id the lattice (universe) id
 * @param num_x the number of lattice cells along x
 * @param num_y the number of lattice cells along y
 * @param origin_x the x-coordinate of the origin
 * @param origin_y the y-coordinate of the origin
 * @param width_x the width of the lattice along x
 * @param width_y the width of the lattice along y
 */
Lattice::Lattice(const int id, const int num_x, int num_y, double origin_x, double origin_y,
		double width_x, double width_y): Universe(id) {
	_num_y = num_y;
	_num_x = num_x;
	_origin.setX(origin_x);
	_origin.setY(origin_y);
	_width_x = width_x;
	_width_y = width_y;
}


/**
 * Lattice destructor
 */
Lattice::~Lattice() {
	for (int i=0; i < _num_x; i++){
		_universes.at(i).clear();
	}
	_universes.clear();
}


/**
 * Add a universe to this lattice
 * @param universe the universe id
 */
void Lattice::addUniverse(const int x, const int y, const int universe) {
	_universes.at(x).at(y) = universe;
	log_printf(INFO, "Added universe with id = %d to lattice with id = %d", universe, _id);
}


/**
 * Get the lattice id
 * @return the lattice id
 */
int Lattice::getId() const {
    return _id;
}


/**
 * Return the number of lattice cells along the x-axis
 * @return the number of lattice cells
 */
int Lattice::getNumX() const {
    return _num_x;
}


/**
 * Return the number of lattice cells along the y-axis
 */
int Lattice::getNumY() const {
    return _num_y;
}


/**
 * Return the origin of the lattice
 * @return the origin of the lattice
 */
Point* Lattice::getOrigin() {
    return &_origin;
}


/**
 * Return a 2D vector array of the universes in the lattice
 * @return 2D vector of universes
 */
std::vector<std::vector<int> > Lattice::getUniverses() const {
    return _universes;
}


/**
 * Return the width of the lattice along the x-axis
 * @return the width of the lattice
 */
double Lattice::getWidthX() const {
    return _width_x;
}


/**
 * Return the width of the lattice along the y-axis
 * @return the width of the lattice
 */
double Lattice::getWidthY() const {
    return _width_y;
}


/**
 * Converts a lattice's attributes to a character array representation
 * @return character array of this lattice's attributes
 */
const char* Lattice::toString() {
	std::stringstream string;

	string << "Material id = " << _id << " num cells along x = "
			<< _num_x << " num cells along y = " << _num_y << " x width = "
			<< _width_x << " y width = " << _width_y;

	string << "\n\tUniverse ids within this lattice:\n\t";
	for (int U = 0; U < (int)_universes.size();  U++) {
		for (int u = 0; u < (int)_universes.at(U).size(); u++)
			string << _universes.at(U).at(u) << "  ";
		string << "\n\t\t";
	}

	string << std::endl;

	return string.str().c_str();
}
