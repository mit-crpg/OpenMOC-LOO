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
#include <string>
#include <math.h>
#include "Parser.h"
#include "Material.h"
#include "Surface.h"
#include "Cell.h"
#include "Universe.h"
#include "Lattice.h"
#include "LocalCoords.h"
#include "Track.h"
#include "log.h"
#include "configurations.h"

class Geometry {
private:
	double _x_min, _y_min, _x_max, _y_max; 		/* the corners */
	int _num_sectors, _num_rings;
	double _sector_offset;
	int _base_universe;
	int _num_FSRs;
	std::map<int, Material*> _materials;
	std::map<int, Surface*> _surfaces;
	std::map<int, Cell*> _cells;
	std::map<int, Universe*> _universes;
	std::map<int, Lattice*> _lattices;

public:
	Geometry(int num_sectors, int num_rings, double sector_offset, 
			Parser* parser);
	virtual ~Geometry();
	void setNumRings(int num_rings);
	void setNumSectors(int num_sectors);
	void setSectorOffset(double sector_offset);
	double getWidth() const;
	double getHeight() const;
	int getNumRings() const;
	int getNumSectors() const;
	double getSectorOffset() const;
	int getNumFSRs() const;

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
	std::string toString();
	void printString();

	void adjustKeys();
	void buildNeighborsLists();
	Cell* findCell(LocalCoords* coords);
	Cell* findNextCell(LocalCoords* coords, double angle);
	void segmentize(Track* track);

	template <class K, class V>
	bool mapContainsKey(std::map<K, V> map, K key);

	void checkUniverse();

};

#endif /* GEOMETRY_H_ */
