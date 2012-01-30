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

class LocalCoords {
private:
	int _cell;
	int _universe;
	int _lattice;
	int _lattice_x;
	int _lattice_y;
	Point _coords;
	LocalCoords* next;
public:
	LocalCoords(double x, double y);
	virtual ~LocalCoords();
    int getCell() const;
    int getUniverse() const;
    int getLattice() const;
    int getLatticeX() const;
    int getLatticeY() const;
    double getX() const;
    double getY() const;
    Point* getPoint();
    LocalCoords *getNext() const;
    void setCell(int cell);
    void setUniverse(int universe);
    void setLattice(int lattice);
    void setLatticeX(int lattice_x);
    void setLatticeY(int lattice_y);
    void setX(double x);
    void setY(double y);
    void setNext(LocalCoords *next);
    std::string toString();
};

#endif /* LOCALCOORDS_H_ */
