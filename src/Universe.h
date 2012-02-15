/*
 * Universe.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "Point.h"
#include "Cell.h"
#include "LocalCoords.h"
#include "log.h"
#include "silo.h"
#include "Surface.h"

class LocalCoords;
class Cell;

enum universeType{
	SIMPLE,
	LATTICE
};

class Universe {
protected:
	static int _n;		/* Counts the number of universes */
	int _uid;		/* monotonically increasing id based on n */
	int _id;
	universeType _type;
	std::map<int, Cell*> _cells;
	Point _origin;
	std::map<int, int> _region_map;
public:
	Universe(const int id);
	virtual ~Universe();
	void addCell(Cell* cell);
	std::map<int, Cell*> getCells() const;
	int getUid() const;
	int getId() const;
	universeType getType();
	int getNumCells() const;
	int getFSR(int cell_id);
	Point* getOrigin();
	void setId(const int id);
	void setType(universeType type);
	void setNumCells(const int num_cells);
	void setOrigin(Point* origin);
	virtual Cell* findCell(LocalCoords* coords,
			       std::map<int, Universe*> universes);
	virtual void generateCSGLists(std::vector<int>* surf_flags, std::vector<double>* surf_coeffs,
			std::vector<int>* oper_flags, std::vector<int>* left_ids, std::vector<int>* right_ids,
			std::vector<int>* zones, Point* point_cur);
	std::string toString();

	virtual int computeFSRMaps();
};

#endif /* UNIVERSE_H_ */

