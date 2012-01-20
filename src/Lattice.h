/*
 * Lattice.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
#include "Point.h"
#include "Universe.h"

class Lattice: public Universe {
private:
	int _id;                       // The universe id
	int _level;                    // The universe level
	int _num_x;						// Number of lattice cells in x direction
	int _num_y;						// Number of lattice cells in y direction
	Point _origin;					// translation origin (global coordinates within core)
	double _width_x;
	double _width_y;
	std::vector< std::vector<int> > _universes;
public:
	Lattice(int id, int num_x, int num_y, double origin_x, double origin_y, double width_x, double width_y);
	virtual ~Lattice();
	void addUniverse(int x, int y, int universe);
    int getId() const;
    int getLevel() const;
    int getNumX() const;
    int getNumY() const;
    Point getOrigin() const;
    std::vector<std::vector<int> > getUniverses() const;
    double getWidthX() const;
    double getWidthY() const;
};

#endif /* LATTICE_H_ */
