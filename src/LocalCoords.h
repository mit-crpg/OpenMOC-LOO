/*
 * LocalCoords.h
 *
 *  Created on: Jan 25, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef LOCALCOORDS_H_
#define LOCALCOORDS_H_

#include <sstream>
#include <string>
#include "Point.h"
#include "Universe.h"

enum coordType {
	UNIV,
	LAT
};

class LocalCoords {
private:
	coordType _type;
	int _universe;
	int _cell;
	int _lattice;
	int _lattice_x;
	int _lattice_y;
	Point _coords;
	LocalCoords* _next;
	LocalCoords* _prev;
public:
	LocalCoords(double x, double y);
	virtual ~LocalCoords();
	coordType getType();
    int getUniverse() const;
    int getCell() const;
    int getLattice() const;
    int getLatticeX() const;
    int getLatticeY() const;
    double getX() const;
    double getY() const;
    Point* getPoint();
    LocalCoords *getNext() const;
    LocalCoords *getPrev() const;
    void setType(coordType type);
    void setUniverse(int universe);
    void setCell(int cell);
    void setLattice(int lattice);
    void setLatticeX(int lattice_x);
    void setLatticeY(int lattice_y);
    void setX(double x);
    void setY(double y);
    void setNext(LocalCoords *next);
    void setPrev(LocalCoords* coords);
    LocalCoords* getLowestLevel();
    void adjustCoords(double delta_x, double delta_y);
    void updateMostLocal(Point* point);
    std::string toString();
};

#endif /* LOCALCOORDS_H_ */
