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
	_type = SIMPLE;
}


/**
 * Destructor
 */
Universe::~Universe() {
	_cells.clear();
}


/**
 * Returns the universe's uid
 * @return the universe's uid
 */
int Universe::getUid() const {
	return _uid;
}


/**
 * Adds a cell to this universe
 * @param cell the cell id
 */
void Universe::addCell(Cell* cell) {
	try {
		_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
		log_printf(INFO, "Added cell with id = %d to universe with id = %d",
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
std::map<int, Cell*> Universe::getCells() const {
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
 * Return the universe type (SIMPLE or LATTICE)
 * @return the universe type
 */
universeType Universe::getType() {
	return _type;
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
 * Sets the universe type to SIMPLE or LATTICE
 * @param type the universe type
 */
void Universe::setType(universeType type) {
	_type = type;
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
 * Finds the cell that a localcoord object is located inside by checking
 * each of this universe's cells. Returns NULL if the localcoord is not
 * in any of the cells
 * @param coords a pointer to the localcoords of interest
 * @param universes a container of all of the universes passed in by geometry
 * @return a pointer the cell where the localcoords is located
 */
Cell* Universe::findCell(LocalCoords* coords,
						std::map<int, Universe*> universes) {

	Cell* return_cell = NULL;
	std::map<int, Cell*>::iterator iter;

	/* Sets the localcoord type to UNIV at this level */
	coords->setType(UNIV);

	/* Loop over all cells in this universe */
	for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
		Cell* cell = iter->second;

		if (cell->cellContains(coords)) {

			/* Set the cell on this level */
			coords->setCell(cell->getId());

			/* MATERIAL type cell - lowest level, terminate search for cell */
			if (cell->getType() == MATERIAL) {
				coords->setCell(cell->getId());
				return_cell = cell;
				return return_cell;
			}

			/* FILL type cell - cell contains a universe at a lower level
			 * Update coords to next level and continue search */
			else if (cell->getType() == FILL) {
				LocalCoords* new_coords = new LocalCoords(coords->getX(),
														coords->getY());
				CellFill* cell_fill = static_cast<CellFill*>(cell);
				int universe_id = cell_fill->getUniverseFill();
				new_coords->setUniverse(universe_id);
				Universe* univ = universes.at(universe_id);
				coords->setCell(cell->getId());

				coords->setNext(new_coords);
				new_coords->setPrev(coords);
				return univ->findCell(new_coords, universes);
			}
		}
	}
	return return_cell;
}


/**
 * Convert the member attributes of this universe to a character array
 * @param a character array reprsenting the universe
 */
std::string Universe::toString() {
	std::stringstream string;

	string << "Universe id = " << _id << ", type = ";
	if (_type == SIMPLE)
		string << "SIMPLE";
	else
		string << "LATTICE";

	string << ", num cells = " << _cells.size() << ", cell ids = ";

	for (int c = 0; c < (int)_cells.size(); c++)
		string << _cells.at(c)->getId() << ", ";

	return string.str();
}
