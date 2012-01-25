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

#include <vector>
#include "log.h"

enum cellType {
	MATERIAL,
	FILL
};

/**
 * Represents a cell
 */
class Cell {
protected:
	int _id;
	cellType _type;
	int _universe;             // universe this cell is in
//	int _parent_cell;          // cell within which this cell resides
	int _num_surfaces;//number of surfaces
	std::vector<int> _surfaces;// + or - depending on side of surface
public:
	Cell(int id, cellType type, int universe, int num_surfaces, 
	     std::vector<int> surfaces);
	virtual ~Cell();
	void addSurface(int surface);
	int getId() const;
	cellType getType() const;
	int getUniverse() const;
//	int getParentCell() const;
	int getNumSurfaces() const;
	std::vector<int> getSurfaces() const;
//	void setParentCell(int parentCell);
	void setUniverse(int universe);
};

/**
 * Represents a cell defined using a material type as a Cell subclass
 */
class CellBasic: public Cell {
private: 
	int _material;             // material filling this cell 	
public:
	CellBasic(int id, int universe, int num_surfaces, 
		  std::vector<int> surfaces, int material);
	int getMaterial() const;
	void setMaterial(int material);
};

/**
 * Represents a cell filled with a universe as a Cell subclass
 */
class CellFill: public Cell {
private:
	int _universe_fill;        // universe filling this cell
public:
	CellFill(int id, int univese, int num_surfaces, 
		 std::vector<int> surfaces, int universe_fill);
	int getUniverseFill() const;
	void setUniverseFill(int universe_Fill);
};

#endif /* CELL_H_ */
