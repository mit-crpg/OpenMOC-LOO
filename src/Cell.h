/*
 * Cell.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef CELL_H_
#define CELL_H_

#include <map>
#include <utility>
#include <sstream>
#include <string>
#include "Surface.h"
#include "Universe.h"
#include "log.h"
#include "Point.h"
#include "LocalCoords.h"

class Universe;

class Surface;
class LocalCoords;

/* Represents cell type */
enum cellType {
    MATERIAL,
    FILL
};

/**
 * Represents a cell inside of a universe
 * @param _n static counter of the number of cells instantiated
 * @param _uid monotonically increasing unique id
 * @param _id cell id from geometry input file
 * @param _type Cell type (MATERIAL or FILL)
 * @param _universe id for the universe this cell is in
 * @param _surfaces map of surface ids (keys) to surface pointers (values)
 */
class Cell {
protected:
    static int _n;
    int _uid;
    int _id;
    cellType _type;
    int _universe;
    std::map<int, Surface*> _surfaces;  /* +/- depending on side of surface */
public:
    Cell();
    Cell(int id, cellType type, int universe, int num_surfaces, 
         int *surfaces);
    virtual ~Cell();
    void addSurface(int surface_id, Surface* surface);
    void setSurfacePointer(Surface* surface);
    int getUid() const;
    int getId() const;
    void setId(int id); 
    cellType getType() const;
    int getUniverse() const;
    int getNumSurfaces() const;
    std::map<int, Surface*> getSurfaces() const;
    void setUniverse(int universe);
    bool cellContains(Point* point);
    bool cellContains(LocalCoords* coords);
    double minSurfaceDist(Point* point, double angle, 
                          Point* min_intersection);
    virtual std::string toString() =0;
    virtual int getNumFSRs() =0;
};


/**
 * Represents a cell filled with a material as a Cell subclass
 * @param _material id for the material filling this cell
 */
class CellBasic: public Cell {
private: 
    int _material;
    int _num_rings;
    int _num_sectors;
public:
    CellBasic(int id, int universe, int num_surfaces, int *surfaces, 
              int material, int num_rings, int num_sectors);
    CellBasic(int id, int universe, int material, 
              int num_rings, int num_sectors);
    CellBasic(int id, int universe, int material);
    int getMaterial() const;
    void addSurface(int surface_id, Surface* surface);
    void adjustKeys(int universe, int material);
    std::string toString();
    int getNumFSRs();
    CellBasic* clone(int new_id, int num_rings, int num_sectors);
    int getNumRings();
    int getNumSectors();
    void setNumSectors(int num);
};


/**
 * Represents a cell filled with a universe as a Cell subclass
 * @param _universe_fill id for the universe filling this cell
 */
class CellFill: public Cell {
private:
    std::pair<int, Universe*> _universe_fill;
public:
    CellFill(int id, int universe, int num_surfaces,
             int *surfaces, int universe_fill);
    int getUniverseFillId() const;
    Universe* getUniverseFill() const;
    void setUniverseFill(int universe_Fill);
    void setUniverseFillPointer(Universe* universe_fill);
    void adjustKeys(int universe, int universe_fill);
    std::string toString();
    int getNumFSRs();
};

#endif /* CELL_H_ */
