/*
 * Universe.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include <vector>
#include <map>
#include <sstream>
#include <string>
#include "Point.h"
#include "Cell.h"
#include "LocalCoords.h"
#include "log.h"

class LocalCoords;

class Cell;

enum universeType{
	SIMPLE,
	LATTICE
};

class Universe {
protected:
	static int _n;				/* Counts the number of universes */
	int _uid;					/* monotonically increasing id based on n */
	int _id;
	universeType _type;
	std::vector<Cell*> _cells;
	Point _origin;
public:
	Universe(const int id);
	virtual ~Universe();
	void addCell(Cell* cell);
    std::vector<Cell*> getCells() const;
    int getUid() const;
    int getId() const;
    universeType getType();
    int getNumCells() const;
    Point* getOrigin();
    void setId(const int id);
    void setType(universeType type);
    void setNumCells(const int num_cells);
    void setOrigin(Point* origin);
    virtual Cell* findCell(LocalCoords* coords, std::map<int, Universe*> universes);
    std::string toString();
};

#endif /* UNIVERSE_H_ */

