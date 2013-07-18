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
#include "LocalCoords.h"
#include "log.h"
#include "silo.h"

class Lattice: public Universe {
private:
    int _num_x, _num_y;
    double _width_x, _width_y;
    std::vector< std::vector< std::pair<int, Universe*> > > _universes;
    std::vector< std::vector< std::pair<int, int> > > _region_map;
    friend class Universe;
	
public:
    Lattice(const int id, const int num_x, const int num_y, 
#ifdef USE_LATTICE_ORIGIN
            const double origin_x, const double origin_y,
#endif
            const double width_x, const double width_y,
            int universes_count, int *universes);
    virtual ~Lattice();
    void setUniversePointer(Universe* universe);
    int getNumX() const;
    int getNumY() const;
    Point* getOrigin();
    std::vector< std::vector< std::pair<int, Universe*>>> 
        getUniverses() const;
    Universe* getUniverse(int lattice_x, int lattice_y) const;
    double getWidthX() const;
    double getWidthY() const;
    int getFSR(int lat_x, int lat_y);
    void adjustKeys();
    bool withinBounds(Point* point);
    Cell* findCell(LocalCoords* coords, std::map<int, Universe*> universes);
    Cell* findNextLatticeCell(LocalCoords* coords, 
                              double angle,
                              std::map<int, Universe*> universes);
    virtual void generateCSGLists(std::vector<int>* surf_flags, std::vector<double>* surf_coeffs,
                                  std::vector<int>* oper_flags, std::vector<int>* left_ids, std::vector<int>* right_ids,
                                  std::vector<int>* zones, Point* point_cur);
    std::string toString();
	
    int virtual computeFSRMaps();
};
#endif /* LATTICE_H_ */
