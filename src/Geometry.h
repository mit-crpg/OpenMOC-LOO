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
#include "LocalCoords.h"
#include "log.h"
#include "configurations.h"

class Geometry {
private:
	double _width, _height;
	int _num_sectors;
	int _num_rings;
	double _sector_offset;
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
    void setNumRings(int num_rings);
    void setNumSectors(int num_sectors);
    void setSectorOffset(double sector_offset);
	double getWidth() const;
	double getHeight() const;
    int getNumRings() const;
    int getNumSectors() const;
    double getSectorOffset() const;

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

	void adjustKeys();
	// Recursively sets Universe levels
//	void buildUniverseLevels(Universe* univ, int parent_cell, int level);
	void buildNeighborsLists();
	bool cellContains(Cell* cell, Point* point);
	bool cellContains(Cell* cell, LocalCoords* coords);

	template <class K, class V>
	bool mapContainsKey(std::map<K, V> map, K key);
};

#endif /* GEOMETRY_H_ */
