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
 * @param universes 
 */
Lattice::Lattice(const int id, const int num_x, int num_y, 
		 double origin_x, double origin_y,
		 double width_x, double width_y, 
		 int universes_count, int *universes): Universe(id) {
	_num_y = num_y;
	_num_x = num_x;
	_origin.setX(origin_x);
	_origin.setY(origin_y);
	_width_x = width_x;
	_width_y = width_y;
	_type = LATTICE;

	Universe* empty_universe_pointer;
	for (int i = 0; i < num_y; i++) {
		_universes.push_back(std::vector< std::pair<int, Universe*> >());
		for (int j = 0; j< num_x; j++){
			_universes.at(i).push_back(std::pair<int, Universe*>
			(universes[i*num_x+j], empty_universe_pointer));
		}
	}
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
void Lattice::setUniversePointer(Universe* universe) {
	/* Check that _surfaces contains this surface id and delete the id
	 *  otherwise
	 * throw an error
	 */
	bool universe_not_found = true;
	int universe_id = universe->getId();

	for (int i = 0; i < _num_y; i++) {
		for (int j = 0; j< _num_x; j++) {
			if (_universes.at(i).at(j).first == universe_id)
				_universes[i][j].second = universe;
			universe_not_found = false;
		}
	}

	if (universe_not_found)
		log_printf(WARNING, "Tried to set the universe pointer for lattice "
				"id = %d for universe id = %d but the lattice does not contain"
				"the universe", _id, universe_id);
	else
		log_printf(INFO, "Set the universe pointer for lattice "
				"id = %d for universe id = %d", _id, universe_id);

	return;
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
std::vector< std::vector< std::pair<int, Universe*> > > Lattice::getUniverses() const {
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
 * Adjusts the ids of the universes inside this lattice to be the uids of each
 * rather than the ids defined by the input file
 */
void Lattice::adjustKeys() {

	std::vector< std::vector< std::pair<int, Universe*> > > adjusted_universes;

	try {
		/* Adjust the indices for each universe to be the universe's uid */
		for (int i = 0; i < _num_y; i++) {
			adjusted_universes.push_back(std::vector< std::pair<int, Universe*> >());

			for (int j = 0; j < _num_x; j++) {
				Universe* universe = _universes.at(i).at(j).second;
				adjusted_universes.at(i).push_back(std::pair<int, Universe*>
												(universe->getUid(), universe));
			}
			_universes.at(i).clear();
		}
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust the keys for lattice id = %s. "
				"Backtrace:\n%s", _id, e.what());
	}

	_universes.clear();
	_universes = adjusted_universes;
}


/**
 * Converts a lattice's attributes to a character array representation
 * @return character array of this lattice's attributes
 */
std::string Lattice::toString() {
	std::stringstream string;

	string << "Lattice id = " << _id << ", num cells along x = "
			<< _num_x << ", num cells along y = " << _num_y << ", x width = "
			<< _width_x << ", y width = " << _width_y;

	string << "\n\t\tUniverse ids within this lattice:\n\t\t";
	for (int i = 0; i < _num_y;  i++) {
		for (int j = 0; j < _num_x; j++)
			string << _universes.at(i).at(j).first << "  ";
		string << "\n\t\t";
	}

	return string.str().c_str();
}
