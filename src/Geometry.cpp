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
		LOG(log_level, "Cannot add a second material with id = %d\n"
				"Exiting program\n", material->getId());
		exit(1);
	}
	else {
		_materials.insert(std::pair<int, Material*>(material->getId(), material));
		LOG(1, "Added material with id = %d to geometry\n", material->getId());
	}
}


/**
 * Return a material from the geometry
 * @param id the material id
 * @return a pointer to the material object
 */
Material* Geometry::getMaterial(int id) {
	if (mapContainsKey(_materials, id))
		return _materials.at(id);
	else {
		LOG(log_level, "Attempted to retrieve material with id = %d which does "
				"not exist\nExiting program\n", id);
		exit(1);
	}
}


/**
 * Add a surface to the geometry
 * @param a pointer to the surface object
 */
void Geometry::addSurface(Surface* surface) {
	if (mapContainsKey(_surfaces, surface->getId())) {
		LOG(log_level, "Cannot add a second surface with id = %d\n"
				"Exiting program\n", surface->getId());
		exit(1);
	}
	else {
		_surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
		LOG(1, "Added surface with id = %d to geometry\n", surface->getId());
	}
}


/**
 * Return a surface from the geometry
 * @param id the surface id
 * @return a pointer to the surface object
 */
Surface* Geometry::getSurface(int id) {
	if (mapContainsKey(_surfaces, id))
		return _surfaces.at(id);
	else {
		LOG(log_level, "Attempted to retrieve surface with id = %d which has"
				" not been declared\nExiting program\n", id);
		exit(1);
	}
}


/**
 * Add a cell to the geometry. Checks if the universe the cell is in already exists;
 * if not, it creates one and adds it to the geometry.
 * @param cell a pointer to the cell object
 */
void Geometry::addCell(Cell* cell) {
	/* If a cell with the same id already exists */
	if (mapContainsKey(_cells, cell->getId())) {
		LOG(log_level, "Cannot add a second cell with id = %d\n"
				"Exiting program\n", cell->getId());
		exit(1);
	}

	/* If the cell's material does not exist */
	else if (cell->getMaterial() != -1E5 && !mapContainsKey(_materials, cell->getMaterial())) {
		LOG(log_level, "Attempted to create cell with material with id = %d, but"
				"material does not exist\nExiting program\n", cell->getMaterial());
		exit(1);
	}

	/* Checks whether the cell's surfaces exist */
	for (int i=0; i < cell->getNumSurfaces(); i++) {
		if (mapContainsKey(_surfaces, abs(cell->getSurfaces().at(i)))) {
			LOG(log_level, "Attempted to create cell with surface id = %d, but "
					"surface does not exist\nExiting Program\n",
					cell->getSurfaces().at(i));
			exit(1);
		}
	}

	_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
	LOG(1, "Added cell with id = %d to geometry\n", cell->getId());

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
	if (mapContainsKey(_cells, id))
		return _cells.at(id);
	else {
		LOG(log_level, "Attempted to retrieve cell with id = %d which has not been "
				"declared\n Exiting program\n", id);
		exit(1);
	}
}


/**
 * Add a universe to the geometry
 * @param universe a pointer to the universe object
 */
void Geometry::addUniverse(Universe* universe) {
	if (mapContainsKey(_universes, universe->getId())) {
		LOG(log_level, "Cannot add a second universe with id = %d\n"
				"Exiting program\n", universe->getId());
		exit(1);
	}
	else {
		_universes.insert(std::pair<int, Universe*>(universe->getId(), universe));
		LOG(1, "Added universe with id = %d to geometry\n", universe->getId());
	}
}


/**
 * Return a universe from the geometry
 * @param the universe id
 * @return a pointer to the universe object
 */
Universe* Geometry::getUniverse(int id) {
	if (mapContainsKey(_universes, id))
		return _universes.at(id);
	else {
		LOG(log_level, "Attempted to retrieve universe with id = %d which has"
				"not been declared\nExiting program\n", id);
		exit(1);
	}
}


/**
 * Add a lattice to the geometry. Adds the lattice to both the lattice and universe containers
 * @param lattice a pointer to the lattice object
 *
 */
void Geometry::addLattice(Lattice* lattice) {
	/* If the lattices container already has a lattice with the same id */
	if (mapContainsKey(_lattices, lattice->getId())) {
		LOG(log_level, "Cannot add a second lattice with id = %d\n"
				"Exiting program\n", lattice->getId());
		exit(1);
	}
	/* If the universes container already has a universe with the same id */
	else if(mapContainsKey(_universes, lattice->getId())) {
		LOG(log_level, "Cannot add a second universe (lattice) with id = %d\n"
				"Exiting program\n", lattice->getId());
		exit(1);
	}
	/* If the lattice contains a universe which does not exist */
	for (int i = 0; i < lattice->getNumX(); i++) {
		for (int j = 0; j < lattice->getNumY(); j++) {
			if (!mapContainsKey(_universes, lattice->getUniverses().at(i).at(j)))
				LOG(log_level, "Attempted to create lattice containing universe"
						"with id = %d, but universe does not exist\nExiting program\n",
						lattice->getUniverses().at(i).at(j));
				exit(1);
		}
	}

	_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
	_universes.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
	LOG(1, "Added lattice with id = %d to geometry\n", lattice->getId());
}


/**
 * Return a lattice from the geometry
 * @param the lattice (universe) id
 * @return a pointer to the lattice object
 */
Lattice* Geometry::getLattice(int id) {
	if (mapContainsKey(_lattices, id))
		return _lattices.at(id);
	else {
		LOG(log_level, "Attempted to retrieve lattice with id = %d which has"
				"not been declared\nExiting program\n", id);
		exit(1);
	}
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
