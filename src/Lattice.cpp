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
#ifdef USE_LATTICE_ORIGIN
	 double origin_x, double origin_y,
#endif
	 double width_x, double width_y,
	 int universes_count, int *universes): Universe(id) {
	_num_y = num_y;
	_num_x = num_x;
#ifdef USE_LATTICE_ORIGIN
	_origin.setX(origin_x);
	_origin.setY(origin_y);
#else
	_origin.setX(-width_x*num_x/2.0);
	_origin.setY(-width_y*num_y/2.0);
#endif
	_width_x = width_x;
	_width_y = width_y;
	_type = LATTICE;

	Universe* empty_universe_pointer;

	/* The parser gives the lattice cells in row major order starting from the
	 * upper left corner. This double loop reorders the lattice cells from the
	 * to start from the lower left corner */
	for (int i = 0; i < num_y; i++) {
		_universes.push_back(std::vector< std::pair<int, Universe*> >());
		for (int j = 0; j< num_x; j++){
			_universes.at(i).push_back(std::pair<int, Universe*>
			(universes[(num_y-1-i)*num_x+j], empty_universe_pointer));
		}
	}

	/* intialize _region_map */
	for (int i = 0; i < num_y; i++) {
		_region_map.push_back(std::vector< std::pair<int, int> >());
		for (int j = 0; j < num_x; j++) {
			_region_map.at(i).push_back(std::pair<int, int>());
		}
	}
}


/**
 * Lattice destructor clears memory for all of its universes
 */
Lattice::~Lattice() {
	for (int i=0; i < _num_x; i++){
		_universes.at(i).clear();
	}
	_universes.clear();
}


/**
 * Sets the pointers to universes in this lattices universe
 * container
 * @param universe the universe id
 */
void Lattice::setUniversePointer(Universe* universe) {
	/* Check that _surfaces contains this surface id and delete the id
	 *  otherwise
	 * throw an error
	 */
	int universe_id = universe->getId();
	bool universe_not_found = true;

	/* Loop over all universes in lattice */
	for (int i = 0; i < _num_y; i++) {
		for (int j = 0; j< _num_x; j++) {
			/* Checks if the universe id matches what the lattice's
			 * universes container expects */
			if (_universes.at(i).at(j).first == universe_id) {
				_universes[i][j].second = universe;
				universe_not_found = false;
			}
		}
	}

	if (universe_not_found)
		log_printf(WARNING, "Tried to set the universe pointer for "
				"lattice id = %d for universe id = %d but the lattice "
				"does not contain the universe", _id, universe_id);
	else
		log_printf(INFO, "Set the universe pointer for lattice "
				"id = %d for universe id = %d", _id, universe_id);

	return;
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
std::vector< std::vector< std::pair<int, Universe*>>> Lattice::getUniverses()
																	const {
    return _universes;
}


/**
 * Returns a universe pointer for a specific lattice cell
 * @param lattice_x the x index to the lattice cell
 * @param lattice_y the y index to the lattice cell
 */
Universe* Lattice::getUniverse(int lattice_x, int lattice_y) const {

	/* Checks that lattice indices are within the bounds of the lattice */
	if (lattice_x > _num_x || lattice_y > _num_y)
		log_printf(ERROR, "Cannot retrieve universe from lattice id = %d: "
				"Index out of bounds: Tried to access cell x = %d, y = %d "
				"but bounds are x = %d, y = %d", _id, lattice_x, lattice_y,
				_num_x, _num_y);

	return _universes.at(lattice_y).at(lattice_x).second;
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


int Lattice::getFSR(int lat_x, int lat_y) {
	/* Check if lattice indices are out of bounds */
	if (lat_x > _num_x || lat_y > _num_y)
		log_printf(ERROR, "Tried to access FSR map of lattice id = %d, but "
				"indices lat_x = %d and lat_y = %d were out of bounds", _id,
				lat_x, lat_y);

	return _region_map[lat_y][lat_x].second;
}


/**
 * Adjusts the ids of the universes inside this lattice to be the uids of each
 * rather than the ids defined by the input file
 */
void Lattice::adjustKeys() {

	/* A new vector of vectors for with universes stored with uids as keys */
	std::vector< std::vector< std::pair<int, Universe*> > > adjusted_universes;

	try {
		/* Adjust the indices for each universe to be the universe's uid */
		for (int i = 0; i < _num_y; i++) {
			adjusted_universes.push_back(std::vector< std::pair<int,
															Universe*> >());

			for (int j = 0; j < _num_x; j++) {
				Universe* universe = _universes.at(i).at(j).second;
				adjusted_universes.at(i).push_back(std::pair<int, Universe*>
											(universe->getUid(), universe));
			}

			/* Delete old universes indexed by ids */
			_universes.at(i).clear();
		}
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust the keys for lattice id = %d. "
				   "Backtrace:\n%s", _id, e.what());
	}

	/* Delete old universes index by ids and set universes container to the
	 * adjusted universes indexed by uids */
	_universes.clear();
	_universes = adjusted_universes;
}


/**
 * Checks if a point is within the bounds of a lattice
 * @param point a pointer to the point of interest
 * @return true if the point is in the bounds, false if not
 */
bool Lattice::withinBounds(Point* point) {

	/* Computes the lattice bounds */
	double bound_x = _num_x/2.0 * _width_x;
	double bound_y = _num_y/2.0 * _width_y;

	double x = point->getX();
	double y = point->getY();

	/* If the point is outside the x bounds */
	if (x > bound_x || x < -1*bound_x)
		return false;
	/* If the point is outside the y boounds */
	else if (y > bound_y || y < -1*bound_y)
		return false;
	/* If the point is within the bounds */
	else
		return true;
}


/**
 * Finds the cell within this lattice that a localcoord is in. This method
 * first find the lattice cell, then searches the universe inside that
 * lattice cell. If localcoord is outside the bounds of the lattice, this
 * method will return NULL
 * @param coords the localcoords of interest
 * @param universes a map of all universes passed in from the geometry
 * @return a pointer to the cell this localcoord is in or NULL
 */
Cell* Lattice::findCell(LocalCoords* coords,
						std::map<int, Universe*> universes) {

	/* Set the localcoord to be a LAT type at this level */
	coords->setType(LAT);

	/* Compute the x and y indices for the lattice cell this coord is in */
	int lat_x = (int)floor((coords->getX() - _origin.getX()) / _width_x);
	int lat_y = (int)floor((coords->getY() - _origin.getY()) / _width_y);

	/* Check if the localcoord is on the lattice boundaries and if so adjust
	 * x or y lattice cell indices i */
	if (fabs(fabs(coords->getX()) - _num_x*_width_x*0.5) <
													ON_LATTICE_CELL_THRESH) {

		if (coords->getX() > 0)
			lat_x = _num_x - 1;
		else
			lat_x = 0;
	}
	if (fabs(fabs(coords->getY()) - _num_y*_width_y*0.5) <
													ON_LATTICE_CELL_THRESH) {
		if (coords->getY() > 0)
			lat_y = _num_y - 1;
		else
			lat_y = 0;
	}

	/* If the indices are outside the bound of the lattice */
	if (lat_x < 0 || lat_x >= _num_x ||
			lat_y < 0 || lat_y >= _num_y) {
		return NULL;
	}

	/* Compute local position of particle in the next level universe */
	double nextX = coords->getX() - (_origin.getX()
					+ (lat_x + 0.5) * _width_x);
	double nextY = coords->getY() - (_origin.getY()
					+ (lat_y + 0.5) * _width_y);

	/* Create a new localcoords object for the next level universe */
	LocalCoords* next_coords;

	if (coords->getNext() == NULL)
		next_coords = new LocalCoords(nextX, nextY);
	else
		next_coords = coords->getNext();

	int universe_id = getUniverse(lat_x, lat_y)->getId();
	Universe* univ = universes.at(universe_id);
	next_coords->setUniverse(universe_id);

	/* Set lattice indices */
	coords->setLattice(_id);
	coords->setLatticeX(lat_x);
	coords->setLatticeY(lat_y);

	coords->setNext(next_coords);
	next_coords->setPrev(coords);

	/* Search the next lowest level universe for the cell */
	return univ->findCell(next_coords, universes);
}


/**
 * Finds the next cell for a localcoords object along a trajectory defined
 * by some angle (in radians from 0 to PI). The method will update the
 * localcoord passed in as an argument to be the one at the boundary of the
 * next cell crossed along the given trajectory. It will do this by
 * recursively building a linked list of localcoords from the localcoord
 * passed in as an argument down to the lowest level cell found. In the
 * process it will set the local coordinates for each localcoord in the
 * linked list for the lattice or universe that it is in. If the
 * localcoord is outside the bounds of the lattice or on the boundaries this
 * method will return NULL; otherwise it will return a pointer to the cell
 * that the localcoords will reach next along its trajectory.
 * @param coords pointer to a localcoords object
 * @param angle the angle of the trajectory
 * @param universes a map of all of the universes passed in by the geometry
 * @return a pointer to a cell if found, NULL if no cell found
 */
Cell* Lattice::findNextLatticeCell(LocalCoords* coords, double angle,
		std::map<int, Universe*> universes) {

	/* Tests the upper, lower, left and right lattice cells adjacent to
	 * the localcoord and uses the one with the shortest distance from
	 * the current location of the localcoord */

	double distance = INFINITY;		/* Initial distance is infinity */
	double d;						/* Current minimum distance */

	/* Properties of the current localcoords */
	double x0 = coords->getX();
	double y0 = coords->getY();
	int lattice_x = coords->getLatticeX();
	int lattice_y = coords->getLatticeY();
	double m = sin(angle) / cos(angle);		/* Slope of trajectory */

	/* Properties of the new location for localcoords */
	double x_curr, y_curr;			/* Current point of minimum distance */
	double x_new = x0;				/* x-coordinate on new lattice cell */
	double y_new = x0; 				/* y-coordinate on new lattice cell */
	int new_lattice_x;	/* New x lattice cell index */
	int new_lattice_y;	/* New y lattice cell index */
	Point test;						/* Test point for computing distance */

	/* Check lower lattice cell Lower lattice cell */
	if (lattice_y >= 0 && angle >= M_PI) {
		y_curr = (lattice_y - _num_y/2.0) * _width_y;
		x_curr = x0 + (y_curr - y0) / m;
		test.setCoords(x_curr, y_curr);

		/* Check if the test point is within the bounds of the lattice */
		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			/* Check if distance to test point is current minimum */
			if (d < distance) {
				distance = d;
				x_new = x_curr;
				y_new = y_curr;
			}
		}
	}

	/* Upper lattice cell */
	if (lattice_y <= _num_y-1 && angle <= M_PI) {
		y_curr = (lattice_y - _num_y/2.0 + 1) * _width_y;
		x_curr = x0 + (y_curr - y0) / m;
		test.setCoords(x_curr, y_curr);

		/* Check if the test point is within the bounds of the lattice */
		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			/* Check if distance to test point is current minimum */
			if (d < distance) {
				distance = d;
				x_new = x_curr;
				y_new = y_curr;
			}
		}
	}

	/* Left lattice cell */
	if (lattice_x >= 0 && (angle >= M_PI/2 && angle <= 3*M_PI/2)) {
		x_curr = (lattice_x - _num_x/2.0) * _width_x;
		y_curr = y0 + m * (x_curr - x0);
		test.setCoords(x_curr, y_curr);

		/* Check if the test point is within the bounds of the lattice */
		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			/* Check if distance to test point is current minimum */
			if (d < distance) {
				distance = d;
				x_new = x_curr;
				y_new = y_curr;
			}
		}
	}

	/* Right lattice cell */
	if (lattice_x <= _num_x-1 && (angle <= M_PI/2 || angle >= 3*M_PI/2)) {
		x_curr = (lattice_x - _num_x/2.0 + 1) * _width_x;
		y_curr = y0 + m * (x_curr - x0);
		test.setCoords(x_curr, y_curr);

		/* Check if the test point is within the bounds of the lattice */
		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			/* Check if distance to test point is current minimum */
			if (d < distance) {
				distance = d;
				x_new = x_curr;
				y_new = y_curr;
			}
		}
	}

	/* If no point was found on the lattice cell, then the localcoords was
	 * already on the boundary of the lattice */
	if (distance == INFINITY)
		return NULL;

	/* Otherwise a point was found inside a new lattice cell */
	else {
		/* Update the localcoords location to the point on the new lattice cell
		 * plus a small bit to ensure that its coordinates are inside cell */
		double delta_x = (x_new - coords->getX()) + cos(angle) * TINY_MOVE;
		double delta_y = (y_new - coords->getY()) + sin(angle) * TINY_MOVE;
		coords->adjustCoords(delta_x, delta_y);

		/* Compute the x and y indices for the new lattice cell */
		new_lattice_x = (int)floor((coords->getX() - _origin.getX()) / _width_x);
		new_lattice_y = (int)floor((coords->getY() - _origin.getY()) / _width_y);

		/* Check if the localcoord is on the lattice boundaries and if so adjust
		 * x or y lattice cell indices i */
		if (fabs(fabs(coords->getX()) - _num_x*_width_x*0.5) <
														ON_LATTICE_CELL_THRESH) {

			if (coords->getX() > 0)
				new_lattice_x = _num_x - 1;
			else
				new_lattice_x = 0;
		}
		if (fabs(fabs(coords->getY()) - _num_y*_width_y*0.5) <
														ON_LATTICE_CELL_THRESH) {
			if (coords->getY() > 0)
				new_lattice_y = _num_y - 1;
			else
				new_lattice_y = 0;
		}

		/* Check if new lattice cell indices are within the bounds, if not, then
		 * new localcoords is now on the boundary of the lattice */
		if (new_lattice_x >= _num_x || new_lattice_x < 0)
			return NULL;
		else if (new_lattice_y >= _num_y || new_lattice_y < 0)
			return NULL;
		/* New localcoords is still within the interior of the lattice */
		else {
			/* Update the localcoords lattice cell indices */
			coords->setLatticeX(new_lattice_x);
			coords->setLatticeY(new_lattice_y);

			/* Move to next lowest level universe */
			coords->prune();
			Universe* univ = _universes.at(new_lattice_y).at(new_lattice_x).second;
			LocalCoords* next_coords;

			/* Compute local position of particle in next level universe */
			double nextX = coords->getX() - (_origin.getX()
							+ (new_lattice_x + 0.5) * _width_x);
			double nextY = coords->getY() - (_origin.getY()
							+ (new_lattice_y + 0.5) * _width_y);

			/* Set the coordinates at the next level localcoord */
			next_coords = new LocalCoords(nextX, nextY);
			next_coords->setPrev(coords);
			coords->setNext(next_coords);

			next_coords->setUniverse(univ->getId());

			/* Search lower level universe */
			return findCell(coords, universes);
		}
	}
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
	for (int i = _num_y-1; i > -1;  i--) {
		for (int j = 0; j < _num_x; j++)
			string << _universes.at(i).at(j).first << "  ";
		string << "\n\t\t";
	}

	return string.str().c_str();
}

int Lattice::computeFSRMaps() {
	/* initialize a counter count */
	int count = 0;
    
	/* loop over universes in the lattice to set the map and update count */
//	for (int i = _num_y -1; i>-1; i--) {
	for (int i = 0; i < _num_y; i++) {
		for (int j = 0; j < _num_x; j++) {
			Universe *u = _universes.at(i).at(j).second;
			_region_map[i][j].second = count;
			count += u->computeFSRMaps();
		}
	}

	return count;
}

void Lattice::generateCSGLists(std::vector<int>* surf_flags, std::vector<double>* surf_coeffs,
		std::vector<int>* oper_flags, std::vector<int>* left_ids, std::vector<int>* right_ids,
		std::vector<int>* zones, Point* current_origin) {

	double global_x, global_y;
	global_x = current_origin->getX();
	global_y = current_origin->getY();

	/* save current cell fill origin */
	Point first_origin;
	first_origin.setCoords(current_origin->getX(), current_origin->getY());

	/* Loop over all lattice cells from left to right, bottom to top */
	for (int i=0; i < _num_y; i++) {

		current_origin->setCoords(current_origin->getX(), global_y + (i * _width_y) - ((float(_num_y)/2.0 - .5) * _width_y));


		for (int j=0; j < _num_x; j++) {

			/*
			 * Get the first cell in this universe. If it is a cell fill, recursively call this function.
			 * If it is a cell basic, create csg zones.
			 */

			Universe* univ;
			univ = getUniverse(j,i);

			std::map<int, Cell*> cells;
			cells = univ->getCells();

			Cell* cell = cells.begin()->second;

			/* update current_origin */
			current_origin->setCoords(global_x + (j * _width_x) - ((float(_num_x)/2.0 - .5) * _width_x),
					current_origin->getY());

			log_printf(DEBUG, "Lattice: set current_origin: (%f,%f)", current_origin->getX(), current_origin->getY());


			/* check cell type */
			if (cell->getType() == MATERIAL){

				/* draw rectangle */

				log_printf(DEBUG, "Lattice: making box (xmin, ymin, xmax, ymax): (%f,%f,%f,%f)",
						current_origin->getX() - _width_x/2.0, current_origin->getX() + _width_x/2.0,
						current_origin->getY() - _width_y/2.0, current_origin->getY() + _width_y/2.0);


				int surf_flags_index = surf_flags->size();

				/* add to surf_flags */
				surf_flags->push_back(DBCSG_LINE_X);
				surf_flags->push_back(DBCSG_LINE_X);
				surf_flags->push_back(DBCSG_LINE_Y);
				surf_flags->push_back(DBCSG_LINE_Y);

				/* add to surf_coeffs */
				surf_coeffs->push_back(current_origin->getX() - _width_x/2.0);
				surf_coeffs->push_back(current_origin->getX() + _width_x/2.0);
				surf_coeffs->push_back(current_origin->getY() - _width_y/2.0);
				surf_coeffs->push_back(current_origin->getY() + _width_y/2.0);

				/* add to oper_flags */
				oper_flags->push_back(DBCSG_OUTER);
				oper_flags->push_back(DBCSG_INNER);
				oper_flags->push_back(DBCSG_OUTER);
				oper_flags->push_back(DBCSG_INNER);
				oper_flags->push_back(DBCSG_INTERSECT);
				oper_flags->push_back(DBCSG_INTERSECT);
				oper_flags->push_back(DBCSG_INTERSECT);

				int left_ids_index = left_ids->size();

				/* add to left_ids */
				left_ids->push_back(surf_flags_index);
				left_ids->push_back(surf_flags_index + 1);
				left_ids->push_back(surf_flags_index + 2);
				left_ids->push_back(surf_flags_index + 3);
				left_ids->push_back(left_ids_index);
				left_ids->push_back(left_ids_index + 2);
				left_ids->push_back(left_ids_index + 4);

				/* add to right_ids */
				right_ids->push_back(-1);
				right_ids->push_back(-1);
				right_ids->push_back(-1);
				right_ids->push_back(-1);
				right_ids->push_back(left_ids_index + 1);
				right_ids->push_back(left_ids_index + 3);
				right_ids->push_back(left_ids_index + 5);

				std::map<int, Surface*> cells_surfaces = cell->getSurfaces();
				std::map<int, Surface*>::iterator iter2;
				double radius;

				/* draw circle */
				for (iter2 = cells_surfaces.begin(); iter2 != cells_surfaces.end(); ++iter2) {

					if (iter2->second->getType() == CIRCLE && iter2->first < 0) {
						radius = static_cast<Circle*>(iter2->second)->getRadius();
						log_printf(DEBUG, "Lattice: making circle (center_x, center_y, radius): (%f,%f,%f)", current_origin->getX(), current_origin->getY(), radius);

						int surf_flags_index = surf_flags->size();

						/* add to surf_flags */
						surf_flags->push_back(DBCSG_CIRCLE_PR);

						/* add to surf_coeffs */
						surf_coeffs->push_back(current_origin->getX());
						surf_coeffs->push_back(current_origin->getY());
						surf_coeffs->push_back(radius);


						/* add to oper_flags */
						oper_flags->push_back(DBCSG_OUTER);
						oper_flags->push_back(DBCSG_INTERSECT);
						oper_flags->push_back(DBCSG_INNER);

						int left_ids_cur = left_ids->size();

						/* add to left_ids */
						left_ids->push_back(surf_flags_index);
						left_ids->push_back(left_ids_cur - 1);
						left_ids->push_back(surf_flags_index);

						/* add to right_ids */
						right_ids->push_back(-1);
						right_ids->push_back(left_ids_cur);
						right_ids->push_back(-1);


						/* add to zones */
						zones->push_back(left_ids_cur + 1);
						zones->push_back(left_ids_cur + 2);

					}
				}
			}

			/* if cell is fill type, recursively call this function */
			else {

				Universe* univ_fill = static_cast<CellFill*>(cell)->getUniverseFill();
//				univ_fill = getUniverse(j,i);

				log_printf(DEBUG, "Lattice: nested universe cell is fill type...");
				univ_fill->generateCSGLists(surf_flags, surf_coeffs, oper_flags,
						left_ids, right_ids, zones, current_origin);
			}
		}
	}

	/* reset current_origin to cell fill origin */
	current_origin->setCoords(first_origin.getX(), first_origin.getY());
	log_printf(DEBUG, "Lattice: set current_origin: (%f,%f)", current_origin->getX(), current_origin->getY());
}


