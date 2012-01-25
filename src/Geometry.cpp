/*
 * Geometry.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Geometry.h"
#include <stdexcept>


/**
 * Geometry constructor
 */
Geometry::Geometry() { }


/**
 * Destructor
 */
Geometry::~Geometry() {
	_materials.clear();
	_surfaces.clear();
	_cells.clear();
	_universes.clear();
	_lattices.clear();
}


/**
 * Sets the total height of the geometry
 * @param height the total height
 */
void Geometry::setHeight(const double height) {
	_height = height;
}


/**
 * Sets the total width of the geometry
 * @param width the total width
 */
void Geometry::setWidth(const double width) {
    _width = width;
}


/**
 * Returns the total height of the geometry
 * @param the toal height of the geometry
 */
double Geometry::getHeight() const {
	return _height;
}


/**
 * Returns the total width of the geometry
 * @param the total width of the geometry
 */
double Geometry::getWidth() const {
    return _width;
}



/**
 * Add a material to the geometry
 * @param material a pointer to a material object
 */
void Geometry::addMaterial(Material* material) {
	if (mapContainsKey(_materials, material->getId())) {
		log_printf(ERROR, "Cannot add a second material with id = %d", material->getId());

	}
	else {
		try {
			_materials.insert(std::pair<int, Material*>(material->getId(), material));
			log_printf(INFO, "Added material with id = %d to geometry\n", material->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add material with id = %d. Backtrace:\n%s",
					material->getId(), e.what());
		}
	}
}


/**
 * Return a material from the geometry
 * @param id the material id
 * @return a pointer to the material object
 */
Material* Geometry::getMaterial(int id) {
	try {
		return _materials.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve material with id = %d which does "
				"not exist. Backtrace:%s\n", id, e.what());
	}
	exit(0);
}


/**
 * Add a surface to the geometry
 * @param a pointer to the surface object
 */
void Geometry::addSurface(Surface* surface) {
	if (mapContainsKey(_surfaces, surface->getId())) {
		log_printf(ERROR, "Cannot add a second surface with id = %d", surface->getId());
	}
	else {
		try {
			_surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
			log_printf(INFO, "Added surface with id = %d to geometry\n", surface->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add surface with id = %d. Backtrace:\n%s",
					surface->getId(), e.what());
		}
	}
}


/**
 * Return a surface from the geometry
 * @param id the surface id
 * @return a pointer to the surface object
 */
Surface* Geometry::getSurface(int id) {
	try {
		return _surfaces.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve surface with id = %d which has"
				" not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}

/**
 * Add a cell to the geometry. Checks if the universe the cell is in already exists;
 * if not, it creates one and adds it to the geometry.
 * @param cell a pointer to the cell object
 */
void Geometry::addCell(Cell* cell) {

	/* If a cell with the same id already exists */
	if (mapContainsKey(_cells, cell->getId()))
		log_printf(ERROR, "Cannot add a second cell with id = %d\n", cell->getId());

	/* If the cell is filled with a material which does not exist */
	else if (cell->getType() == MATERIAL &&
			!mapContainsKey(_materials, static_cast<CellBasic*>(cell)->getMaterial())) {
		log_printf(ERROR, "Attempted to create cell with material with id = %d, but "
			"material does not exist", static_cast<CellBasic*>(cell)->getMaterial());
	}

	/* If the cell is filled with a universe which doesn't exist yet, create it */
	else if (cell->getType() == FILL && !mapContainsKey(_universes,
			static_cast<CellFill*>(cell)->getUniverse())) {

		Universe* univ = new Universe(cell->getUniverse());
		addUniverse(univ);
	}


	/* Checks whether the cell's surfaces exist */
	for (int i=0; i < cell->getNumSurfaces(); i++) {
		if (!mapContainsKey(_surfaces, abs(cell->getSurfaces().at(i))))
			log_printf(ERROR, "Attempted to create cell with surface id = %d, but "
					"surface does not exist", cell->getSurfaces().at(i));
	}

	/* Insert the cell into the geometry's cell container */
	try {
		_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
		log_printf(INFO, "Added cell with id = %d to geometry\n", cell->getId());
	}
	catch (std::exception &e) {
			log_printf(ERROR, "Unable to add cell with id = %d. Backtrace:\n%s",
					cell->getId(), e.what());
	}

	/* Checks if the universe the cell in exists and if not, creates a new universe */
	if (!mapContainsKey(_universes, cell->getUniverse())) {
		try {
			Universe* univ = new Universe(cell->getUniverse());
			addUniverse(univ);
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to create a new universe with id = %d and add "
					"it to the geometry. Backtrace:\n%s", cell->getUniverse(), e.what());
		}
	}

}


/**
 * Return a cell from the geometry
 * @param id the cell's id
 * @return a pointer to the cell object
 */
Cell* Geometry::getCell(int id) {
	try {
		return _cells.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve cell with id = %d which has not been "
				"declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a universe to the geometry
 * @param universe a pointer to the universe object
 */
void Geometry::addUniverse(Universe* universe) {
	if (mapContainsKey(_universes, universe->getId()))
		log_printf(ERROR, "Cannot add a second universe with id = %d", universe->getId());
	else {
		try {
			_universes.insert(std::pair<int, Universe*>(universe->getId(), universe));
			log_printf(INFO, "Added universe with id = %d to geometry\n", universe->getId());
		}
		catch (std::exception &e) {
				log_printf(ERROR, "Unable to add universe with id = %d. Backtrace:\n%s",
						universe->getId(), e.what());
		}
	}
}


/**
 * Return a universe from the geometry
 * @param the universe id
 * @return a pointer to the universe object
 */
Universe* Geometry::getUniverse(int id) {
	try {
		return _universes.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve universe with id = %d which has not been "
				"declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a lattice to the geometry. Adds the lattice to both the lattice and universe containers
 * @param lattice a pointer to the lattice object
 *
 */
void Geometry::addLattice(Lattice* lattice) {
	/* If the lattices container already has a lattice with the same id */
	if (mapContainsKey(_lattices, lattice->getId()))
		log_printf(ERROR, "Cannot add a second lattice with id = %d", lattice->getId());

	/* If the universes container already has a universe with the same id */
	else if(mapContainsKey(_universes, lattice->getId()))
		log_printf(ERROR, "Cannot add a second universe (lattice) with id = %d", lattice->getId());

	/* If the lattice contains a universe which does not exist */
	for (int i = 0; i < lattice->getNumX(); i++) {
		for (int j = 0; j < lattice->getNumY(); j++) {
			if (!mapContainsKey(_universes, lattice->getUniverses().at(i).at(j)))
				log_printf(ERROR, "Attempted to create lattice containing universe"
						"with id = %d, but universe does not exist",
						lattice->getUniverses().at(i).at(j));
		}
	}

	try {
		_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
		log_printf(INFO, "Added lattice with id = %d to geometry\n", lattice->getId());
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add lattice with id = %d. Backtrace:\n%s",
				lattice->getId(), e.what());
	}

	/* Add the lattice to the universes container as well */
	addUniverse(lattice);
}


/**
 * Return a lattice from the geometry
 * @param the lattice (universe) id
 * @return a pointer to the lattice object
 */
Lattice* Geometry::getLattice(int id) {
	try {
		return _lattices.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve lattice with id = %d which has"
				"not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Converts this geometry's attributes to a character array
 * @param a character array of this geometry's attributes
 */
const char* Geometry::toString() {
	std::stringstream string;
	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;


	string << "Geometry: width = " << _width << ", height = " << _height <<
			", base universe id = " << _base_universe;

	string << "\n\tMaterial ids = ";
	for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
		string << iter1->second << ", ";

	string << "\n\tSurface ids = ";
	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
		string << iter2->second << ", ";

	string << "\n\tCell ids = ";
	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
		string << iter3->second << ", ";

	string << "\n\tUniverse ids = ";
	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4)
		string << iter4->second << ", ";

	string << "\n\tLattice ids = ";
	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5)
		string << iter5->second << ", ";

	string << std::endl;

	return string.str().c_str();
}


// Adjusts the keys for surfaces, cells, universes, and lattices to uids
void Geometry::adjustKeys() {

	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;

	std::map<int, Material*> adjusted_materials;
	std::map<int, Surface*> adjusted_surfaces;
	std::map<int, Cell*> adjusted_cells;
	std::map<int, Universe*> adjusted_universes;
	std::map<int, Lattice*> adjusted_lattices;

	int uid;


	/**************************************************************************
	 * Ajust the indices for attributes of all cell, universe and lattice
	 * objects in the geometry
	 *************************************************************************/

	/* Adjust the container of surface ids inside each cell to hold the surfaces' uids */
	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3) {

		Cell* cell = iter3->second;
		int universe = _universes.at(cell->getUniverse())->getUid();

		/* MATERIAL type cells */
		if (cell->getType() == MATERIAL) {
			CellBasic* cell_basic = static_cast<CellBasic*>(cell);
			int material = _materials.at(cell_basic->getMaterial())->getUid();
			cell_basic->adjustKeys(universe, material, _surfaces);
		}

		/* FILL type cells */
		else {
			CellFill* cell_fill = static_cast<CellFill*>(cell);
			int universe_fill = _universes.at(cell_fill->getUniverseFill())->getUid();
			cell_fill->adjustKeys(universe, universe_fill, _surfaces);
		}
	}

	/* Adjust the container of cell ids inside each cell to hold the cells' uids */
	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4) {
		Universe* universe = iter4->second;
		universe->adjustKeys(_cells);
	}

	/* Adjust the container of universe ids inside each lattice to hold the universes' uids */
	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
		Lattice* lattice = iter5->second;
		lattice->adjustKeys(_universes);
	}


	/**************************************************************************
	 * Ajust the indices of the containers of geometry objects which are
	 * attributes of this geometry class
	 *************************************************************************/

	/* Adjust material indices to be uids */
	try {
		for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1) {
			uid = iter1->second->getUid();
			Material* material = iter1->second;
			adjusted_materials.insert(std::pair<int, Material*>(uid, material));
		}
		_materials.clear();
		_materials = adjusted_materials;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust material' keys. Backtrace:\n%s", e.what());
	}


	/* Adjust surfaces indices to be uids */
	try {
		for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2) {
			uid = iter2->second->getUid();
			Surface* surface = iter2->second;
			adjusted_surfaces.insert(std::pair<int, Surface*>(uid, surface));
		}
		_surfaces.clear();
		_surfaces = adjusted_surfaces;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust surface' keys. Backtrace:\n%s", e.what());
	}

	/* Adjust cells indices to be uids */
	try {
		for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3) {
			uid = iter3->second->getUid();
			Cell* cell = iter3->second;
			adjusted_cells.insert(std::pair<int, Cell*>(uid, cell));
		}
		_cells.clear();
		_cells = adjusted_cells;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust cell' keys Backtrace:\n%s", e.what());
	}

	/* Adjust universes indices to be uids */
	try {
		for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4) {
			uid = iter4->second->getUid();
			Universe* universe = iter4->second;
			adjusted_universes.insert(std::pair<int, Universe*>(uid, universe));
		}
		_universes.clear();
		_universes = adjusted_universes;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust universes' keys Backtrace:\n%s", e.what());
	}

	/* Adjust lattices indices to be uids */
	try {
		for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
			uid = iter5->second->getUid();
			Lattice* lattice = iter5->second;
			adjusted_universes.insert(std::pair<int, Lattice*>(uid, lattice));
		}
		_lattices.clear();
		_lattices = adjusted_lattices;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust lattices' keys Backtrace:\n%s", e.what());
	}

	return;
}



/**
 * Function to determine whether a key already exists in a templated map container
 * @param map the map container
 * @param key the key to check
 * @return true if already in map, false otherwise
 */
template <class K, class V>
bool Geometry::mapContainsKey(std::map<K, V> map, K key) {
	/* Try to access the element at the key */
	try { map.at(key); }

	/* If an exception is thrown, element does not exist */
	catch (std::exception& exc) { return false; }

	/* If no exception is thrown, element does exist */
	return true;
}
