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
#include <utility>
#include <string>
#include "Point.h"
#include "Universe.h"
#include "log.h"

class Lattice: public Universe {
private:
	int _num_x, _num_y;
	double _width_x, _width_y;
	std::vector< std::vector< std::pair<int, Universe*> > > _universes;
	friend class Universe;
public:
	Lattice(const int id, const int num_x, const int num_y, 
		const double origin_x, const double origin_y,
		const double width_x, const double width_y,
		int universes_count, int *universes);
	virtual ~Lattice();
	void setUniversePointer(Universe* universe);
    int getId() const;
    int getLevel() const;
    int getNumX() const;
    int getNumY() const;
    Point* getOrigin();
    std::vector< std::vector< std::pair<int, Universe*> > > getUniverses() const;
    double getWidthX() const;
    double getWidthY() const;
    void adjustKeys();
    std::string toString();
};

#endif /* LATTICE_H_ */
