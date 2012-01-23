/*
 * Geometry.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Geometry.h"

/**
 * Geometry constructor
 */
Geometry::Geometry() {
}

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
 * Add a material to the geometry
 * @param material a pointer to a material object
 */
void Geometry::addMaterial(Material* material) {
	if (mapContainsKey(_materials, material->getId())) {
		log_printf(ERROR, "Cannot add a second material with id = %d", material->getId());

	}
	else {
		_materials.insert(std::pair<int, Material*>(material->getId(), material));
		log_printf(INFO, "Added material with id = %d to geometry\n", material->getId());
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
		log_printf(ERROR, "Attempted to retrieve material with id = %d which does not exist", id);
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
		_surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
		log_printf(INFO, "Added surface with id = %d to geometry\n", surface->getId());
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
				" not been declared", id);
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

	/* If the cell's material does not exist */
	else if (cell->getMaterial() != -1E5 && !mapContainsKey(_materials, cell->getMaterial()))
		log_printf(ERROR, "Attempted to create cell with material with id = %d, but "
				"material does not exist", cell->getMaterial());

	/* Checks whether the cell's surfaces exist */
	for (int i=0; i < cell->getNumSurfaces(); i++) {
		if (mapContainsKey(_surfaces, abs(cell->getSurfaces().at(i))))
			log_printf(ERROR, "Attempted to create cell with surface id = %d, but "
					"surface does not exist", cell->getSurfaces().at(i));
	}

	_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
	log_printf(INFO, "Added cell with id = %d to geometry\n", cell->getId());

	/* Checks if the universe the cell in exists and if not, creates a new universe */
	if (!mapContainsKey(_universes, cell->getUniverse())) {
		Universe* univ = new Universe(cell->getUniverse());
		addUniverse(univ);
	}
}


/**
 * Return a cell from the geometry
 * @param id the cell's id
 * @return a pointer to the cell object
 */
Cell* Geometry::getCell(int id) {
	try {
		log_printf(ERROR, "Attempted to retrieve surface with id = %d which has"
				" not been declared", id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve cell with id = %d which has not been declared", id);
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
		_universes.insert(std::pair<int, Universe*>(universe->getId(), universe));
		log_printf(INFO, "Added universe with id = %d to geometry\n", universe->getId());
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
		log_printf(ERROR, "Attempted to retrieve universe with id = %d which has not been declared", id);
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

	_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
	_universes.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
	log_printf(INFO, "Added lattice with id = %d to geometry\n", lattice->getId());
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
				"not been declared", id);
	}
	exit(0);
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
