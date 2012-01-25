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


class Cell {
private:
	int _id;
	cellType _type;
	int _parent_cell;          // cell within which this cell resides
	int _universe;             // universe this cell is in
	int _universe_fill;        // universe filling this cell
	int _material;             // material filling this cell (0 for base universe)
	int _num_surfaces;
	std::vector<int> _surfaces;     // positive or negative ids depending on side of surface
public:
	Cell(int id, cellType type, int num_surfaces);
	virtual ~Cell();
    void addSurface(int surface);
    int getId() const;
    int getMaterial() const;
    int getNumSurfaces() const;
    int getParentCell() const;
    std::vector<int> getSurfaces() const;
    cellType getType() const;
    int getUniverse() const;
    int getUniverseFill() const;
    void setParentCell(int parentCell);
    void setMaterial(int material);
    void setUniverse(int universe);
    void setUniverseFill(int universeFill);
};

#endif /* CELL_H_ */
