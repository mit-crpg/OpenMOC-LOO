/*
 * Universe.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include <vector>
#include <sstream>
#include <string>
#include "Point.h"
#include "Cell.h"
#include "log.h"

class Universe {
protected:
	static int _n;				/* Counts the number of universes */
	int _uid;					/* monotonically increasing id based on n */
	int _id;
	std::vector<Cell*> _cells;
	Point _origin;
public:
	Universe(const int id);
	virtual ~Universe();
	void addCell(Cell* cell);
    std::vector<Cell*> getCells() const;
    int getUid() const;
    int getId() const;
    int getNumCells() const;
    Point* getOrigin();
    void setId(const int id);
    void setNumCells(const int num_cells);
    void setOrigin(Point* origin);
    std::string toString();
};

#endif /* UNIVERSE_H_ */

