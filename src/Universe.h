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
#include "Point.h"
#include "Cell.h"
#include "log.h"

class Universe {
protected:
	static int _n;				/* Counts the number of universes */
	int _uid;					/* monotonically increasing id based on n */
	int _id;
	int _num_cells;
	std::vector<int> _cells;
	Point _origin;
public:
	Universe(const int id);
	virtual ~Universe();
	void addCell(const int cell);
    std::vector<int> getCells() const;
    int getUid() const;
    int getId() const;
    int getNumCells() const;
    Point* getOrigin();
    void setCells(std::vector<int> cells);
    void setId(const int id);
    void setNumCells(const int numCells);
    void setOrigin(Point* origin);
    void adjustKeys(std::map<int, Cell*> cells);
    const char* toString();
};

#endif /* UNIVERSE_H_ */

