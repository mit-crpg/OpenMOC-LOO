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
		std::cout << "Two materials with the same id = " << material->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else
		_materials.insert(std::pair<int, Material*>(material->getId(), material));
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
		std::cout << "No material with id = " << id << " has been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
}


/**
 * Add a surface to the geometry
 * @param a pointer to the surface object
 */
void Geometry::addSurface(Surface* surface) {
	if (mapContainsKey(_surfaces, surface->getId())) {
		std::cout << "Two surfaces with the same id = " <<  surface->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else
		_surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
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
		std::cout << "No surface with id = " << id << " has been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
}


/**
 * Add a cell to the geometry. Checks if the universe the cell is in already exists;
 * if not, it creates one and adds it to the geometry.
 * @param cell a pointer to the cell object
 */
void Geometry::addCell(Cell* cell) {
	if (mapContainsKey(_cells, cell->getId())) {
		std::cout << "Two cells with the same id = " <<  cell->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else {
		_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));

		// Checks if the universe the cell in exists and if not, creates a new
		// universe
		if (!mapContainsKey(_universes, cell->getUniverse())) {
			Universe* univ = new Universe(cell->getUniverse());
			addUniverse(univ);
		}
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
		std::cout << "No cell with id = " << id << " has been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
}


/**
 * Add a universe to the geometry
 * @param universe a pointer to the universe object
 */
void Geometry::addUniverse(Universe* universe) {
	if (mapContainsKey(_universes, universe->getId())) {
		std::cout << "Two universes with the same id = " <<  universe->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else
		_universes.insert(std::pair<int, Universe*>(universe->getId(), universe));
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
		std::cout << "No universe with id = " << id << " has been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
}


/**
 * Add a lattice to the geometry. Adds the lattice to both the lattice and universe containers
 * @param lattice a pointer to the lattice object
 *
 */
void Geometry::addLattice(Lattice* lattice) {
	if (mapContainsKey(_lattices, lattice->getId())) {
		std::cout << "Two lattices with the same id = " <<  lattice->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else if(mapContainsKey(_universes, lattice->getId())) {
		std::cout << "Two universes with the same id = " <<  lattice->getId() << " have been declared.";
		std::cout << " Exiting program. " << std::endl;
		exit(1);
	}
	else {
		_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
		_universes.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
	}
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
		std::cout << "No lattice with id = " << id << " has been declared.";
		std::cout << " Exiting program. " << std::endl;
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
	// Try to access the element at the key
	try { map.at(key); }

	// If an exception is thrown, element does not exist
	catch (std::exception& exc) { return false; }

	// If no exception is thrown, element does exist
	return true;
}
