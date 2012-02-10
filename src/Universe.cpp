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


int Universe::getFSR(int cell_id) {

	if (_cells.find(cell_id) == _cells.end())
		log_printf(ERROR, "Tried to find FSR id for cell with id = %d was found"
				" in universe with id = %d but no cell exists", cell_id, _id);

	return _region_map.at(cell_id);
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

				LocalCoords* next_coords;

				if (coords->getNext() == NULL)
					next_coords = new LocalCoords(coords->getX(),
													coords->getY());
				else
					next_coords = coords->getNext();

				CellFill* cell_fill = static_cast<CellFill*>(cell);
				int universe_id = cell_fill->getUniverseFillId();
				next_coords->setUniverse(universe_id);
				Universe* univ = universes.at(universe_id);
				coords->setCell(cell->getId());

				coords->setNext(next_coords);
				next_coords->setPrev(coords);
				return univ->findCell(next_coords, universes);
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
	std::map<int, Cell*>::iterator iter;

	string << "Universe id = " << _id << ", type = ";
	if (_type == SIMPLE)
		string << "SIMPLE";
	else
		string << "LATTICE";

	string << ", num cells = " << _cells.size() << ", cell ids = ";

	for (iter = _cells.begin(); iter != _cells.end(); ++iter)
		string << iter->first << ", ";

	return string.str();
}

/*
 * Compute the FSR Maps for this universe
 */
int Universe::computeFSRMaps() {
	/* initialize a counter count */
	std::map<int, Cell*>::iterator iter;
	int count = 0;
    
	/* loop over cells in the universe to set the map and update count */
	for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
		_region_map.insert(std::pair<int, int>(iter->first, count));
		count += iter->second->getNumFSRs();
	}

	return count;
}

void Universe::generateCSGLists(std::vector<int>* surf_flags, std::vector<double>* surf_coeffs,
		std::vector<int>* oper_flags, std::vector<int>* left_ids, std::vector<int>* right_ids,
		std::vector<int>* zones, Point* point, Point* point_cur){

	/* Get the starting index into the vectors passed in...this will be the current length of the vector and may come in handy for when you create surfaces and regions and zones */

	std::map<int, Cell*>::iterator iter;

	/* Loop over all of this universes cells */
	for (iter = _cells.begin(); iter != _cells.end(); ++iter){
		Cell* cell = iter->second;

		log_printf(DEBUG, "Inside Universe CSGLists method for id = %d, cell id = %d", _id, cell->getId());

		if (cell->getType() == MATERIAL) {
			/* Do everything you would do to add surfaces to the vectors and
			 *  make them intersect into regions and create a new zone.
			 *   This will assume that the most recent set of surfaces are
			 *   those boundary planes surrounding the pin cell, and that the
			 *    most recent region created is the box of surfaces surrounding
			 *   the pin cell (see notes below for Lattice.cpp)
			 */

			std::map<int, Surface*> cells_surfaces = cell->getSurfaces();
			std::map<int, Surface*>::iterator iter2;
			double radius;

			for (iter2 = cells_surfaces.begin(); iter2 != cells_surfaces.end(); ++iter2) {

				log_printf(DEBUG, "Checking surface id = %d", iter2->first);
				if (iter2->second->getType() == CIRCLE && iter2->first < 0) {
					radius = static_cast<Circle*>(iter2->second)->getRadius();

					// add to surf_flags
					log_printf(DEBUG, "surf flags size: %d", surf_flags->size());
					surf_flags->push_back(DBCSG_CIRCLE_PR);
					int surf_index = surf_flags->size();
					log_printf(DEBUG, "Circle csg id: %d", DBCSG_CIRCLE_PR);
					log_printf(DEBUG, "box csg id: %d", DBCSG_BOX_XYXY);

					// add to surf_coeffs
					surf_coeffs->push_back(point_cur->getX());
					surf_coeffs->push_back(point_cur->getY());
					surf_coeffs->push_back(radius);

					log_printf(DEBUG, "Checking new surf_coeffs size = %d", surf_coeffs->size());

					// add to oper_flags
					oper_flags->push_back(DBCSG_OUTER);
					oper_flags->push_back(DBCSG_INTERSECT);
					oper_flags->push_back(DBCSG_INNER);

					int left_ids_cur = left_ids->size();
					log_printf(DEBUG, "left_ids_cur: %d", left_ids->size());

					// add to left_ids
					left_ids->push_back(surf_index);
					left_ids->push_back(left_ids_cur - 1);
					left_ids->push_back(surf_index);

					// add to right_ids
					right_ids->push_back(-1);
					right_ids->push_back(left_ids_cur);
					right_ids->push_back(-1);


					// add to zones
					zones->push_back(left_ids_cur + 1);
					zones->push_back(left_ids_cur + 2);

				}
			}

		}

		/* If the Universe contains a FILL type Cell, then recursively go into
		 *  the universe filling this cell
		 */
		else if (cell->getType() == FILL) {

			Universe* univ = static_cast<CellFill*>(cell)->getUniverseFill();
			univ->generateCSGLists(surf_flags, surf_coeffs, oper_flags,
					left_ids, right_ids, zones, point, point_cur);
		}
	}
}


