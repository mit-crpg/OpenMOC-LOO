/*
 * Geometry.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <map>
#include <utility>
#include <sstream>
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Universe.h"
#include "Lattice.h"
#include "log.h"


class Geometry {
private:
	double _width, _height;
	int _base_universe;
	std::map<int, Material*> _materials;
	std::map<int, Surface*> _surfaces;
	std::map<int, Cell*> _cells;
	std::map<int, Universe*> _universes;
	std::map<int, Lattice*> _lattices;               // Needed? Or do we put lattices inside of universes?
public:
	Geometry();
	virtual ~Geometry();
	void setWidth(const double width);
	void setHeight(const double height);
	double getWidth() const;
	double getHeight() const;
	void addMaterial(Material* material);
	Material* getMaterial(int id);
	void addSurface(Surface* surface);
	Surface* getSurface(int id);
	void addCell(Cell *cell);
	Cell* getCell(int id);
	void addUniverse(Universe* universe);
	Universe* getUniverse(int id);
	void addLattice(Lattice* lattice);
	Lattice* getLattice(int id);
	const char* toString();

	// Adjusts the keys for surfaces, cells, universes, and lattices to uids
	void adjustKeys();
	// Recursively sets Universe levels
	void buildUniverseLevels(Universe* univ, int parent_cell, int level);
	// Build list of neighboring cells in positive/negative direction from each surface
	void buildNeighborsList();

	template <class K, class V>
	bool mapContainsKey(std::map<K, V> map, K key);
};

#endif /* GEOMETRY_H_ */
