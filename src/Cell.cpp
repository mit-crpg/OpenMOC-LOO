/*
 * Cell.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Cell.h"

/* _n keeps track of the number cells of instantiated */
int Cell::_n = 0;


/**
 * Default Cell constructor
 */
Cell::Cell() { }

/**
 * Cell constructor
 * @param id the cell id
 * @param type the type of cell
 * @param universe the universe this cell is in
 * @param parent_cell the cell within which this cell resides
 * @param num_surfaces the number of surfaces in this cell
 * @param surfaces the surface id
 */
Cell::Cell(int id, cellType type, int universe, int num_surfaces,
           int *surfaces) {

    _uid = _n;
    _id = id;
    _type = type;
    _universe = universe;
    _n++;
	
    /* This empty surface pointer is just a null value for the _surfaces
     * map. The Geometry will register the actual surface pointer when
     * the cell is added to the geometry */
    Surface* empty_surface_pointer;
    for (int i = 0; i < num_surfaces; i++)
        _surfaces.insert(std::pair<int, Surface*>(surfaces[i],
                                                  empty_surface_pointer));
}


/**
 * Destructor frees all surfaces making up cell
 */
Cell::~Cell() {
    _surfaces.clear();
}

/**
 * Insert the a surface into this cell's container
 * @param surface_id the surface id (positiv/negative for surface side)
 * @param surface surface pointer
 */
void Cell::addSurface(int surface_id, Surface* surface) {
    _surfaces.insert(std::pair<int, Surface*>(surface_id, surface));
}



/**
 * Registers a surface pointer with the Cell's surfaces map. This method is
 * used by the geometry whenever a new cell is added to it
 * @param surface the surface pointer
 */
void Cell::setSurfacePointer(Surface* surface) {

    /* if _surfaces does not contain this surface id throw an error */
    if (_surfaces.find(surface->getId()) == _surfaces.end() &&
        _surfaces.find(-surface->getId()) == _surfaces.end())

        log_printf(WARNING, "Unable to set surface pointer for cell id = %d "
                   "for surface id = %d since cell does not contain this surface",
                   _id, surface->getId());

    try{
        /* If the cell contains the positive side of the surface */
        if (_surfaces.find(surface->getId()) != _surfaces.end())
            _surfaces[surface->getId()] = surface;

        /* If the cell contains the negative side of the surface */
        else
            _surfaces[-1*surface->getId()] = surface;

        log_printf(INFO, "Set the surface pointer for cell id = %d for "
                   "surface id = %d", _id, surface->getId());

        log_printf(INFO, "Set the surface pointer for cell id = %d for "
                   "surface id = %d", _id, surface->getId());

        if (surface->getType() == XPLANE){
            log_printf(INFO, "Surface %i has x = %f", 
                       surface->getId(), surface->getXMin());
        }
        else if (surface->getType() == YPLANE){
            log_printf(INFO, "Surface %i has y = %f", 
                       surface->getId(), surface->getYMin());
        }

    }
    catch (std::exception &e) {
        log_printf(ERROR, 
                   "Unable to add surface with id = %d to cell with id = %d. "
                   "Backtrace:\n%s", surface->getId(), _id, e.what());
    }
}


/**
 * Return the cell's uid
 * @return the cell's uid
 */
int Cell::getUid() const {
    return _uid;
}



/**
 * Return the cell's id
 * @return the cell's id
 */
int Cell::getId() const {
    return _id;
}

void Cell::setId(int id) {
    _id = id;
}


/**
 * Return the cell type (FILL or MATERIAL)	Cell();
 * @return the cell type
 */
cellType Cell::getType() const {
    return _type;
}


/**
 * Return the universe that this cell is in
 * @return the universe id
 */
int Cell::getUniverse() const {
    return _universe;
}


/**
 * Return the number of surfaces in the cell
 * @return the number of surfaces
 */
int Cell::getNumSurfaces() const {
    return this->_surfaces.size();
}



/**
 * Return the vector of surfaces in the cell
 * @return vector of surface ids
 */
std::map<int,Surface*> Cell::getSurfaces() const {
    return _surfaces;
}


/**
 * Set the universe that this cell is inside of
 * @param the universe's id
 */
void Cell::setUniverse(int universe) {
    _universe = universe;
}


/*
 * Determines whether a point is contained inside a cell. Queries each surface
 * inside the cell to determine if the particle is on the same side of the
 * surface. This particle is only inside the cell if it is on the same side of
 * every surface in the cell.
 * @param point a pointer to a point
 */
bool Cell::cellContains(Point* point) {

    /* Loop over all surfaces inside the cell */
    std::map<int, Surface*>::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        /* If the surface evaluated point is not the same sign as the surface
         * or within a threshold, return false */
        if (iter->second->evaluate(point) * iter->first < -ON_SURFACE_THRESH){
            return false;
        }

    }

    return true;
}


/*
 * Determines whether a point is contained inside a ce	void clone(int new_id);
 * ll. Queries each surface
 * inside the cell to determine if the particle is on the same side of the
 * surface. This particle is only inside the cell if it is on the same side of
 * every surface in the cell.
 * @param point a pointer to a localcoord
 */
bool Cell::cellContains(LocalCoords* coords) {
    return this->cellContains(coords->getPoint());
}


/*	Cell();
 *
 * Computes the minimum distance to a surface from a point with a given
 * trajectory at a certain angle. If the trajectory will not intersect
 * any of the surfaces in the cell, returns INFINITY
 * @param point the point of interest
 * @param angle the angle of the trajectory (in radians from 0 to 2*PI)
 * @param min_intersection a pointer to the intersection point that is found
 */
double Cell::minSurfaceDist(Point* point, double angle,
                            Point* min_intersection) {

    double min_dist = INFINITY;	void clone(int new_id);

    double d;
    Point intersection;

    std::map<int, Surface*>::iterator iter;

    /* Loop over all of the cell's surfaces */
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        /* Find the minimum distance from this surface to this point */
        d = iter->second->getMinDistance(point, angle, &intersection);

        /* If the distance to cell is less than current min distance, update */
        if (d < min_dist) {
            min_dist = d;
            min_intersection->setX(intersection.getX());
            min_intersection->setY(intersection.getY());
        }
    }

    return min_dist;
}


/**
 *  CellBasic constructor
 *  @param id the cell id
 *  @param universe the id of the universe this cell is in
 *  @param num_surfaces the number of surfaces in this cell
 *  @param surfaces array of surface ids in this cell
 *  @param material id of the material filling this cell
 */

CellBasic::CellBasic(int id, int universe, int num_surfaces, int* surfaces, 
                     int material, int num_rings, int num_sectors):
    Cell(id, MATERIAL, universe, num_surfaces, surfaces) {
    _material = material;
    _num_rings = num_rings;
    _num_sectors = num_sectors;
}


/**
 * Constructor for the CellBasic's clone method
 * @param id the cell id
 * @param universe the id of the universe this cell is in
 * @param num_surfaces the number of surfaces inside this cell
 */
CellBasic::CellBasic(int id, int universe, int material, 
                     int num_rings, int num_sectors){
    _id = id;
    _universe = universe;
    _material = material;
    _num_rings = num_rings;
    _num_sectors = num_sectors;
}

CellBasic::CellBasic(int id, int universe, int material){ 
    _id = id;
    _universe = universe;
    _material = material;
}
/**
 * Return the material in the cell	int getUniverseFillId() const;
 * @return the material's id
 */
int CellBasic::getMaterial() const {
    return _material;
}

/**
 * Insert the a surface into this cell's container
 * @param surface_id the surface id (positiv/negative for surface side)
 * @param surface surface pointer
 */
void CellBasic::addSurface(int surface_id, Surface* surface) {
    _surfaces.insert(std::pair<int, Surface*>(surface_id, surface));
}


/**
 * Adjusts the cell's id of the universe this cell is in, the id of the
 * universe that is filling this cell, and the id of the surfaces inside
 * this cell to be the uids of each rather than the ids defined by the
 * input file
 * @param universe uid of the universe this cell is in
 * @param material uid of the material this cell contains
 */
void CellBasic::adjustKeys(int universe, int material) {

    _universe = universe;
    _material = material;

    /* New surfaces map uses uid instead of id as the keys */
    std::map<int, Surface*> adjusted_surfaces;

    try {
        /* Adjust surface ids to be the positive/negative surface uids */
        std::map<int, Surface*>::iterator iter;
        for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

            /* Positive side of surface, use positive uid */
            if (iter->first > 0)
                adjusted_surfaces.insert(std::pair<int, Surface*>
                                         (iter->second->getUid(), iter->second));

            /* Negative side of surface, use negative uid */
            else
                adjusted_surfaces.insert(std::pair<int, Surface*>
                                         (iter->second->getUid()*-1, iter->second));

        }

        /* Reset the surfaces container to be the one with new keys*/
        _surfaces.clear();
        _surfaces = adjusted_surfaces;
    }
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to adjust the surface ids to uids for cell "
                   "id = %d. Backtrace:\n%s", _id, e.what());
    }
}


/**
 * Convert this cell's attributes to a string format
 * @return a character array of this cell's attributes
 */
std::string CellBasic::toString() {

    std::stringstream string;

    string << "Cell id = " << _id 
           << ", type = MATERIAL, material id = " << _material 
           << ", universe = " << _universe  
           << ", num_surfaces = " << getNumSurfaces() 
           << ", num of rings = " << _num_rings 
           << ", num of sectors = " << _num_sectors;

    std::map<int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


CellBasic* CellBasic::clone(int new_id, int num_rings, int num_sectors) {

    CellBasic* new_cell = new CellBasic(new_id, _universe, _material, 
                                        num_rings, num_sectors);

    /* Loop over all of this cell's surfaces and add them to the clone */
    std::map<int, Surface*>::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        new_cell->addSurface(iter->first, iter->second);
    return new_cell;
}

/**
 * Return the number of rings in the cell
 * @return the number of rings
 */
int CellBasic::getNumRings() {
    return this->_num_rings;
}

/**
 * Return the number of sectors in the cell
 * @return the number of sectors
 */
int CellBasic::getNumSectors() {
    return this->_num_sectors;
}

/**
 * set the cell's number of sectors 
 */
void CellBasic::setNumSectors(int num) {
    _num_sectors = num;
}

/**
 *  CellFill constructor
 *  @param id the cell id
 *  @param universe the id of the universe this cell is in
 *  @param num_surfaces the number of surfaces in this cell
 *  @param surfaces array of surface ids in this cell
 *  @param universe_fill the universe used to fill this cell
 */
CellFill::CellFill(int id, int universe, int num_surfaces, 
		   int *surfaces, int universe_fill):
    Cell(id, FILL, universe, num_surfaces, surfaces) {

    _universe_fill.first = universe_fill;
}

/**
 * Return the universe filling this cell
 * @return the universe's id
 */
int CellFill::getUniverseFillId() const {
    return _universe_fill.first;
}

Universe* CellFill::getUniverseFill() const {
    return _universe_fill.second;
}


/**
 * Set the universe filling this cell
 * @param the universe's id
 */
void CellFill::setUniverseFill(int universe_fill) {
    _universe_fill.first = universe_fill;
}


void CellFill::setUniverseFillPointer(Universe* universe) {
    _universe_fill.second = universe;
}

/**
 * Adjusts the cell's id of the universe this cell is in, the id of the
 * universe that is filling this cell, and the id of the surfaces inside
 * this cell to be the uids of each rather than the ids defined by the
 * input file
 * @param universe uid of the universe this cell is in
 * @param material uid of the material this cell contains
 *
 */
void CellFill::adjustKeys(int universe, int universe_fill) {

    _universe = universe;
    _universe_fill.first = universe_fill;

    /* New surfaces map uses uid instead of id as the keys */
    std::map<int, Surface*> adjusted_surfaces;

    try {
        /* Adjust surface ids to be the positive/negative surface uids */
        std::map<int, Surface*>::iterator iter;
        for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

            /* Positive side of surface, use positive uid */
            if (iter->first > 0)
                adjusted_surfaces.insert(std::pair<int, Surface*>
                                         (iter->second->getUid(), iter->second));

            /* Negative side of surface, use negative uid */
            else
                adjusted_surfaces.insert(std::pair<int, Surface*>
                                         (iter->second->getUid()*-1, iter->second));
        }

        /* Reset the surfaces container to be the one with new keys*/
        _surfaces.clear();
        _surfaces = adjusted_surfaces;
    }
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to adjust the surface ids to uids for cell "
                   "id = %d. Backtrace:\n%s", _id, e.what());
    }

}


/**
 * Convert this cell's attributes to a string format
 * @return a character array of this cell's attributes
 */
std::string CellFill::toString() {

    std::stringstream string;

    string << "Cell id = " << _id << ", type = FILL, universe_fill = " <<
        _universe_fill.first << ", universe = " << _universe << ", num_surfaces = "
           << getNumSurfaces();

    std::map<int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}

int CellBasic::getNumFSRs() {
    return 1;
}

int CellFill::getNumFSRs() {
    Universe *univ = getUniverseFill();
    return univ->computeFSRMaps();
}
