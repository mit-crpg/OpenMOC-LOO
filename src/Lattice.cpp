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
	 * to start from th lower left corner */
	for (int i = 0; i < num_y; i++) {
		_universes.push_back(std::vector< std::pair<int, Universe*> >());
		for (int j = 0; j< num_x; j++){
			_universes.at(i).push_back(std::pair<int, Universe*>
			(universes[(num_y-1-i)*num_x+j], empty_universe_pointer));
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


#ifdef USE_LATTICE_ORIGIN
/**
 * Return the origin of the lattice
 * @return the origin of the lattice
 */
Point* Lattice::getOrigin() {
    return &_origin;
}
#endif

/**
 * Return a 2D vector array of the universes in the lattice
 * @return 2D vector of universes
 */
std::vector< std::vector< std::pair<int, Universe*> > > Lattice::getUniverses() const {
    return _universes;
}


Universe* Lattice::getUniverse(int lattice_x, int lattice_y) const {
	if (lattice_x > _num_x || lattice_y > _num_y)
		log_printf(ERROR, "Cannot retrieve universe from lattice id = %d: Index"
				" out of bounds: Tried to access cell x = %d, y = %d but bounds "
				"are x = %d, y = %d", _id, lattice_x, lattice_y, _num_x, _num_y);

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
 * @return the width of the lattice    Cell* findNextLatticeCell();
 *
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

bool Lattice::withinBounds(Point* point) {

	double bound_x = _num_x/2.0 * _width_x;
	double bound_y = _num_y/2.0 * _width_y;
	double x = point->getX();
	double y = point->getY();

	if (x > bound_x || x < -1*bound_x)
		return false;
	else if (y > bound_y || y < -1*bound_y)
		return false;
	else
		return true;
}


Cell* Lattice::findCell(LocalCoords* coords, std::map<int, Universe*> universes) {

	coords->setType(LAT);

	/* Compute the x and y indices for the lattice cell this coord is in */
	int lat_x = floor(coords->getX() - _origin.getX()) / _width_x;
	int lat_y = floor(coords->getY() - _origin.getY()) / _width_y;

//	int lat_x = ceil((coords->getX() - _origin.getX()) / _width_x);
//	int lat_y = ceil((coords->getY() - _origin.getY()) / _width_y);


//	if (coords->getX() < ON_LATTICE_CELL_THRESH && coords->getX() > -ON_LATTICE_CELL_THRESH) {
//		lat_x = 0;
//		log_printf(DEBUG, "Resetting lat_x since it is within the THRESH");
//	}
//	if (coords->getY() < ON_LATTICE_CELL_THRESH && coords->getY() > -ON_LATTICE_CELL_THRESH) {
//		lat_y = 0;
//		log_printf(DEBUG, "Resetting lat_x since it is within the THRESH");
//	}

	log_printf(DEBUG, "coords.getY() = %1.10f, _origin.getY() = %1.10f, floor(.....) = %f", coords->getY(), _origin.getY(), floor(coords->getY() - _origin.getY())/ _width_y);
	log_printf(DEBUG, "Inside Lattice findCell. x = %1.10f, y = %1.10f, lat_x = %d, lat_y = %d, _width_y = %f", coords->getX(), coords->getY(), lat_x, lat_y, _width_y);

	if (fabs(fabs(coords->getX()) - _num_x*_width_x*0.5) < ON_LATTICE_CELL_THRESH) {
		if (coords->getX() > 0)
			lat_x = _num_x - 1;
		else
			lat_x = 0;

		log_printf(DEBUG, "ON_LATTICE_CELL_THRESH, updated lat_x = %d", lat_x);
	}
	if (fabs(fabs(coords->getY()) - _num_y*_width_y*0.5) < ON_LATTICE_CELL_THRESH) {
		if (coords->getY() > 0)
			lat_y = _num_y - 1;
		else
			lat_y = 0;

		log_printf(DEBUG, "ON_LATTICE_CELL_THRESH, updated lat_y = %d", lat_y);
	}

	/* If the indices are outside the bound of the lattice */
	if (lat_x < 0 || lat_x >= _num_x ||
			lat_y < 0 || lat_y >= _num_y) {

		log_printf(WARNING, "The lattice cell indices are out of "
				"bounds (x = %d, y = %d) for the following "
				"LocalCoords and and Lattice objects:\n%s\n%s",
				lat_x, lat_y, coords->toString().c_str(), toString().c_str());

		return NULL;
	}

	/* Compute local position of particle in the next level universe */
	double nextX = coords->getX() - (_origin.getX()
					+ (lat_x + 0.5) * _width_x);
	double nextY = coords->getY() - (_origin.getY()
					+ (lat_y + 0.5) * _width_y);

	LocalCoords* newCoords = new LocalCoords(nextX, nextY);
	int universe_id = getUniverse(lat_x, lat_y)->getId();
	Universe* univ = universes.at(universe_id);
	newCoords->setUniverse(universe_id);

	/* Set lattice indices */
//	coords->setLattice(_uid);
	coords->setLattice(_id);
	coords->setLatticeX(lat_x);
	coords->setLatticeY(lat_y);

	coords->setNext(newCoords);
	newCoords->setPrev(coords);

	return univ->findCell(newCoords, universes);
}


Cell* Lattice::findNextLatticeCell(LocalCoords* coords, double angle,
		std::map<int, Universe*> universes) {

	double distance = INFINITY;
	double d;
	double x0 = coords->getX();
	double y0 = coords->getY();
	double x1, y1;
	double x_new, y_new;
	double m = sin(angle) / cos(angle);
	int lattice_x = coords->getLatticeX();
	int lattice_y = coords->getLatticeY();
	int new_lattice_x = lattice_x;
	int new_lattice_y = lattice_y;

	log_printf(DEBUG, "Inside latticeFindNextCell. lat_x = %d, lat_y = %d, x0 = %f, y0 = %f", lattice_x, lattice_y, x0, y0);

	Point test;

	/* Lower lattice cell */
	if (lattice_y >= 0 && angle >= M_PI) {
//		y1 = (lattice_y - 1) * _width_y;
//		y1 = (lattice_y - 2) * _width_y;
		y1 = (lattice_y - _num_y/2.0) * _width_y;
		x1 = x0 + (y1 - y0) / m;
		test.setCoords(x1, y1);

		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			if (d < distance) {
				log_printf(DEBUG, "Moving to lower lattice cell. dist = %f, x1 = %f, y1 = %f", d, x1, y1);
				distance = d;
				new_lattice_x = lattice_x;
				new_lattice_y = lattice_y - 1;
				x_new = x1;
				y_new = y1;
			}
		}
	}

	/* Upper lattice cell */
	if (lattice_y <= _num_y-1 && angle <= M_PI) {
//		y1 = lattice_y * _width_y;
//		y1 = (lattice_y - 1) * _width_y;
		y1 = (lattice_y - _num_y/2.0 + 1) * _width_y;
		x1 = x0 + (y1 - y0) / m;
		test.setCoords(x1, y1);

		log_printf(DEBUG, "Testing upper lattice cell. x1 = %f, y1 = %f", x1, y1);

		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			if (d < distance) {
				log_printf(DEBUG, "Moving to upper lattice cell. dist = %f, x1 = %f, y1 = %f", d, x1, y1);
				distance = d;
				new_lattice_x = lattice_x;
				new_lattice_y = lattice_y + 1;
				x_new = x1;
				y_new = y1;
			}
		}
	}

	/* Left lattice cell */
	if (lattice_x >= 0 && (angle >= M_PI/2 && angle <= 3*M_PI/2)) {
//		x1 = (lattice_x - 1) * _width_x;
//		x1 = (lattice_x - 2) * _width_x;
		x1 = (lattice_x - _num_x/2.0) * _width_x;
		y1 = y0 + m * (x1 - x0);
		test.setCoords(x1, y1);

		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			if (d < distance) {
				log_printf(DEBUG, "Moving to left lattice cell. dist = %f, x1 = %f, y1 = %f", d, x1, y1);
				distance = d;
				new_lattice_x = lattice_x -1 ;
				new_lattice_y = lattice_y;
				x_new = x1;
				y_new = y1;
			}
		}
	}

	/* Right lattice cell */
	if (lattice_x <= _num_x-1 && (angle <= M_PI/2 || angle >= 3*M_PI/2)) {
//		x1 = (lattice_x) * _width_x;
//		x1 = (lattice_x - 1) * _width_x;
		x1 = (lattice_x - _num_x/2.0 + 1) * _width_x;
		y1 = y0 + m * (x1 - x0);
		test.setCoords(x1, y1);

		log_printf(DEBUG, "Testing right lattice cell. x1 = %f, y1 = %f", x1, y1);

		if (withinBounds(&test)) {
			d = test.distance(coords->getPoint());

			if (d < distance) {
				log_printf(DEBUG, "Moving to right lattice cell. dist = %f, x1 = %f, y1 = %f", d, x1, y1);
				distance = d;
				new_lattice_x = lattice_x + 1;
				new_lattice_y = lattice_y;
				x_new = x1;
				y_new = y1;
			}
		}
	}

	if (distance == INFINITY) {
		log_printf(DEBUG, "Distance to next lattice cell is INFINITY");
		return NULL;
	}

	else {
		double delta_x = (x_new - coords->getX()) + cos(angle) * TINY_MOVE;
		double delta_y = (y_new - coords->getY()) + sin(angle) * TINY_MOVE;
//		double delta_x = cos(angle) * (distance + TINY_MOVE);
//		double delta_y = sin(angle) * (distance + TINY_MOVE);
		coords->adjustCoords(delta_x, delta_y);

		if (new_lattice_x >= _num_x || new_lattice_x < 0) {
			log_printf(DEBUG, "Returning NULL from findNextLatticeCell method. new_lattice_x = %d, localcoords: %s", new_lattice_x, coords->toString().c_str());
			return NULL;
		}
		else if (new_lattice_y >= _num_y || new_lattice_y < 0) {
			log_printf(DEBUG, "Returning NULL from findNextLatticeCell method. new_lattice_y = %d, localcoords: %s", new_lattice_y, coords->toString().c_str());
			return NULL;
		}
		else {
			coords->setLatticeX(new_lattice_x);
			coords->setLatticeY(new_lattice_y);
			Universe* univ = _universes.at(new_lattice_y).at(new_lattice_x).second;
			LocalCoords* new_coords;

			if (coords->getNext() != NULL)
				new_coords = coords->getNext();
			else {
				/* Compute local position of particle in the next level universe */
				double nextX = coords->getX() - (_origin.getX()
								+ (new_lattice_x + 0.5) * _width_x);
				double nextY = coords->getY() - (_origin.getY()
								+ (new_lattice_y + 0.5) * _width_y);

				new_coords = new LocalCoords(nextX, nextY);
				new_coords->setPrev(coords);
				coords->setNext(new_coords);
				log_printf(DEBUG, "Moving coordinates in lower universe");
			}

			log_printf(DEBUG, "new coords: %s", new_coords->toString().c_str());
			new_coords->setUniverse(univ->getId());

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
