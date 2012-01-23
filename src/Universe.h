/*
 * Universe.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include <vector>
#include "Point.h"
#include "log.h"

class Universe {
private:
	int _id;
	int _level;
	int _num_cells;
	std::vector<int> _cells;
	Point _origin;
public:
	Universe(int id);
	virtual ~Universe();
	void addCell(int cell);
    std::vector<int> getCells() const;
    int getId() const;
    int getLevel() const;
    int getNumCells() const;
    Point* getOrigin();
    void setCells(std::vector<int> cells);
    void setId(int id);
    void setLevel(int level);
    void setNumCells(int numCells);
    void setOrigin(Point* origin);
};

#endif /* UNIVERSE_H_ */

