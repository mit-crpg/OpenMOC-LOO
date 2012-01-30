/*
 * Cell.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef CELL_H_
#define CELL_H_

#include <map>
#include <utility>
#include <sstream>
#include <string>
#include "Surface.h"
#include "log.h"
#include "Point.h"
#include "LocalCoords.h"

class Surface;
class LocalCoords;

enum cellType {
	MATERIAL,
	FILL
};

/**
 * Represents a cell
 */
class Cell {
protected:
	static int _n;				/* Counts the number of cells */
	int _uid;					/* monotonically increasing id based on n */
	int _id;
	cellType _type;
	int _universe;             	/* universe id this cell is in */
//	std::vector<int> _surfaces;	/* + or - depending on side of surface */
	std::map<int, Surface*> _surfaces;
public:
	Cell(int id, cellType type, int universe, int num_surfaces, 
	     int *surfaces);
	virtual ~Cell();
//	void addSurface(int surface);
	void setSurfacePointer(Surface* surface);
	int getUid() const;
	int getId() const;
	cellType getType() const;
	int getUniverse() const;
	int getNumSurfaces() const;
	std::map<int, Surface*> getSurfaces() const;
	void setUniverse(int universe);
	bool cellContains(Point* point);
	bool cellContains(LocalCoords* coords);
	virtual std::string toString() =0;
};


/**
 * Represents a cell defined using a material type as a Cell subclass
 */
class CellBasic: public Cell {
private: 
	int _material;             // material filling this cell 	
public:
	CellBasic(int id, int universe, int num_surfaces, 
		  int *, int material);
	int getMaterial() const;
	void setMaterial(int material);
	void adjustKeys(int universe, int material);
	std::string toString();
};


/**
 * Represents a cell filled with a universe as a Cell subclass
 */
class CellFill: public Cell {
private:
	int _universe_fill;        // universe filling this cell
public:
	CellFill(int id, int universe, int num_surfaces,
		 int *surfaces, int universe_fill);
	int getUniverseFill() const;
	void setUniverseFill(int universe_Fill);
	void adjustKeys(int universe, int universe_fill);
	std::string toString();
};

#endif /* CELL_H_ */
