/*
 * Geometry.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Geometry.h"
#include <stdexcept>


/**
 * Geometry constructor
 * @param num_sectors the number of angular sectors per cell
 * @param num_rings the number of rings per cell
 * @param sector_offset the angular offset for computing angular sectors
 */
Geometry::Geometry(int num_sectors, int num_rings, double sector_offset) {

	_num_sectors = num_sectors;
	_num_rings = num_rings;
	_sector_offset = sector_offset;

	/* Initializing the corners to be infinite  */
	_x_min = 1.0/0.0;
	_y_min = 1.0/0.0;
	_x_max = -1.0/0.0;
	_y_max = -1.0/0.0;

 }


/**
 * Destructor
 */
Geometry::~Geometry() {
	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;

//	for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
//		delete iter1->second;
	_materials.clear();

//	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
//		delete iter2->second;
	_surfaces.clear();

//	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
//		delete iter3->second;
	_cells.clear();

//	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4)
//		delete iter4->second;
	_universes.clear();

//	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5)
//		delete iter5->second;
	_lattices.clear();
}

#if 0
/**
 * Sets the total height of the geometry
 * @param height the total height
 */
void Geometry::setHeight(const double height) {
	_height = height;
}


/**
 * Sets the total width of the geometry
 * @param width the total width
 */
void Geometry::setWidth(const double width) {
    _width = width;
}
#endif


/* Set the number of ring divisions used for making flat source regions
 * @param num_rings the number of rings
 */
void Geometry::setNumRings(int num_rings) {
    _num_rings = num_rings;
}


/**
 * Set the number of angular sectors used for making flat source regions
 * @param num_sectors the number of sectors
 */
void Geometry::setNumSectors(int num_sectors) {
    _num_sectors = num_sectors;
}


/**
 * Sets the angular offset (from the positive x-axis) by which angular
 * divisions are computed. Takes the modulo of the input argument with
 * 360 degrees
 * @param angular_offset angular offset of sectors in degrees
 */
void Geometry::setSectorOffset(double sector_offset) {
    _sector_offset = fmod(sector_offset, 360);
}


/**
 * Returns the total height of the geometry
 * @param the toal height of the geometry
 */
double Geometry::getHeight() const {
	return _y_max - _y_min;
}


/**
 * Returns the total width of the geometry
 * @param the total width of the geometry
 */
double Geometry::getWidth() const {
    return _x_max - _x_min;
}


/**
 * Returns the number of rings used in subdividing flat source regions
 * @return the number of rings
 */
int Geometry::getNumRings() const{
    return _num_rings;
}


/**
 * Returns the number of angular sectors used in subdividing flat source regions
 * @param the number of angular sectors
 */
int Geometry::getNumSectors() const {
    return _num_sectors;
}


/**
 * Returns the angular offset for the angular divisions used in subdividing
 * flat source regions
 * @return the angular offset in degrees
 */
double Geometry::getSectorOffset() const {
    return _sector_offset;
}


/**
 * Add a material to the geometry
 * @param material a pointer to a material object
 */
void Geometry::addMaterial(Material* material) {
	if (_materials.find(material->getId()) != _materials.end())
		log_printf(ERROR, "Cannot add a second material with id = %d", material->getId());

	else {
		try {
			_materials.insert(std::pair<int, Material*>(material->getId(), material));
			log_printf(INFO, "Added material with id = %d to geometry\n", material->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add material with id = %d. Backtrace:\n%s",
					material->getId(), e.what());
		}
	}
}


/**
 * Return a material from the geometry
 * @param id the material id
 * @return a pointer to the material object
 */
Material* Geometry::getMaterial(int id) {
	try {
		return _materials.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve material with id = %d which does "
				"not exist. Backtrace:%s\n", id, e.what());
	}
	exit(0);
}


/**
 * Add a surface to the geometry
 * @param a pointer to the surface object
 */
void Geometry::addSurface(Surface* surface) {
	if (_surfaces.find(surface->getId()) != _surfaces.end())
		log_printf(ERROR, "Cannot add a second surface with id = %d", surface->getId());

	else {
		try {
			_surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
			log_printf(INFO, "Added surface with id = %d to geometry\n", surface->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add surface with id = %d. Backtrace:\n%s",
					surface->getId(), e.what());
		}
	}

	switch (surface->getBoundary()) {
	case REFLECTIVE:
		if (surface->getXMin() <= _x_min)
			_x_min = surface->getXMin();
		if (surface->getXMax() >= _x_max)
			_x_max = surface->getXMax();
		if (surface->getYMin() <= _y_min)
			_y_min = surface->getYMin();
		if (surface->getYMax() >= _y_max)
			_y_max = surface->getYMax();
		break;
	case BOUNDARY_NONE:
		break;
	}
	
}


/**
 * Return a surface from the geometry
 * @param id the surface id
 * @return a pointer to the surface object
 */
Surface* Geometry::getSurface(int id) {
	try {
		return _surfaces.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve surface with id = %d which has"
				" not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a cell to the geometry. Checks if the universe the cell is in already exists;
 * if not, it creates one and adds it to the geometry.
 * @param cell a pointer to the cell object
 */
void Geometry::addCell(Cell* cell) {

	/* If a cell with the same id already exists */
	if (_cells.find(cell->getId()) != _cells.end())
		log_printf(ERROR, "Cannot add a second cell with id = %d\n", cell->getId());

	/* If the cell is filled with a material which does not exist */
	else if (cell->getType() == MATERIAL &&
			_materials.find(static_cast<CellBasic*>(cell)->getMaterial()) ==
														_materials.end()) {

		log_printf(ERROR, "Attempted to add cell with material with id = %d,"
				" but material does not exist",
				static_cast<CellBasic*>(cell)->getMaterial());
	}

	/* If the cell is filled with a universe which doesn't exist, create it */
	else if (cell->getType() == FILL &&
			_universes.find(static_cast<CellFill*>(cell)->getUniverse()) ==
														_universes.end()) {

		Universe* univ = new Universe(cell->getUniverse());
		addUniverse(univ);
	}


	/* Set the pointers for each of the surfaces inside the cell and also
	 * checks whether the cell's surfaces exist */
	std::map<int, Surface*> cells_surfaces = cell->getSurfaces();
	std::map<int, Surface*>::iterator iter;

	/* Loop over all surfaces in the cell */
	for (iter = cells_surfaces.begin(); iter != cells_surfaces.end(); ++iter) {
		int surface_id = abs(iter->first);

		/* The surface does not exist */
		if (_surfaces.find(surface_id) == _surfaces.end())
			log_printf(ERROR, "Attempted to add cell with surface id = %d, "
					"but surface does not exist\n", iter->first);

		/* The surface does exist, so set the surface pointer in the cell */
		else
			cell->setSurfacePointer(_surfaces.at(surface_id));
	}

	/* Insert the cell into the geometry's cell container */
	try {
		_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
		log_printf(INFO, "Added cell with id = %d to geometry\n", cell->getId());
	}
	catch (std::exception &e) {
			log_printf(ERROR, "Unable to add cell with id = %d. Backtrace:\n%s",
					cell->getId(), e.what());
	}

	/* Checks if the universe the cell in exists; if not, creates new universe */
	if (_universes.find(cell->getUniverse()) == _universes.end()) {
		try {
			Universe* univ = new Universe(cell->getUniverse());
			addUniverse(univ);
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to create a new universe with id = %d and"
					" add it to the geometry. Backtrace:\n%s",
					cell->getUniverse(), e.what());
		}
	}
	_universes.at(cell->getUniverse())->addCell(cell);

	return;
}


/**
 * Return a cell from the geometry
 * @param id the cell's id
 * @return a pointer to the cell object
 */
Cell* Geometry::getCell(int id) {
	try {
		return _cells.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve cell with id = %d which has "
				"not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a universe to the geometry
 * @param universe a pointer to the universe object
 */
void Geometry::addUniverse(Universe* universe) {
	if (_universes.find(universe->getId()) != _universes.end())
		log_printf(ERROR, "Cannot add a second universe with id = %d",
				universe->getId());
	else {
		try {
			_universes.insert(std::pair<int, Universe*>(universe->getId(),
														universe));
			log_printf(INFO, "Added universe with id = %d to geometry\n",
					universe->getId());
		}
		catch (std::exception &e) {
				log_printf(ERROR, "Unable to add universe with id = %d. "
						"Backtrace:\n%s", universe->getId(), e.what());
		}
	}
}


/**
 * Return a universe from the geometry
 * @param the universe id
 * @return a pointer to the universe object
 */
Universe* Geometry::getUniverse(int id) {
	try {
		return _universes.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve universe with id = %d which "
				"has not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a lattice to the geometry. Adds the lattice to both the lattice and
 * universe containers
 * @param lattice a pointer to the lattice object
 *
 */
void Geometry::addLattice(Lattice* lattice) {
	/* If the lattices container already has a lattice with the same id */
	if (_lattices.find(lattice->getId()) != _lattices.end())
		log_printf(ERROR, "Cannot add a second lattice with id = %d",
				lattice->getId());

	/* If the universes container already has a universe with the same id */
	else if (_universes.find(lattice->getId()) != _universes.end())
		log_printf(ERROR, "Cannot add a second universe (lattice) with "
				"id = %d", lattice->getId());

	/* Sets the universe pointers for the lattice and checks if the lattice
	 * contains a universe which does not exist */
	for (int i = 0; i < lattice->getNumY(); i++) {
		for (int j = 0; j < lattice->getNumX(); j++) {
			int universe_id = lattice->getUniverses().at(i).at(j).first;

			/* If the universe does not exist */
			if (_universes.find(universe_id) == _universes.end())
				log_printf(ERROR, "Attempted to create lattice containing "
						"universe with id = %d, but universe does not exist",
						lattice->getUniverses().at(i).at(j));

			/* Set the universe pointer */
			else
				lattice->setUniversePointer(_universes.at(universe_id));
		}
	}

	/* Add the lattice to the geometry's lattices container */
	try {
		_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
		log_printf(INFO, "Added lattice with id = %d to geometry\n",
						lattice->getId());
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add lattice with id = %d. Backtrace:\n%s",
				lattice->getId(), e.what());
	}

	/* Add the lattice to the universes container as well */
	addUniverse(lattice);
}


/**
 * Return a lattice from the geometry
 * @param the lattice (universe) id
 * @return a pointer to the lattice object
 */
Lattice* Geometry::getLattice(int id) {
	try {
		return _lattices.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve lattice with id = %d which has"
				"not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Converts this geometry's attributes to a character array
 * @param a character array of this geometry's attributes
 */
std::string Geometry::toString() {

	std::stringstream string;
	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;


	string << "Geometry: width = " << getWidth() << ", height = " << getHeight()
			<< ", base universe id = " << _base_universe << ", Bounding Box: (("
			<< _x_min << ", " << _y_min << "), (" << _x_max << ", " << _y_max
			<< ")";

	string << "\n\tMaterials:\n\t\t";
	for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
		string << iter1->second->toString() << "\t\t";

	string << "\n\tSurfaces:\n\t\t";
	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
		string << iter2->second->toString() << "\t\t";

	string << "\n\tCells:\n\t\t";
	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
		string << iter3->second->toString() << "\t\t";

	string << "\n\tUniverses:\n\t\t";
	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4)
		string << iter4->second->toString() << "\t\t";

	string << "\n\tLattices:\n\t\t";
	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
		string << iter5->second->toString()  << "\t\t";
	}

	string << "\n";
	return string.str();
}


// Adjusts the keys for surfaces, cells, universes, and lattices to uids
void Geometry::adjustKeys() {

	log_printf(NORMAL, "Adjusting the keys for the geometry...\n");

	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;

	std::map<int, Material*> adjusted_materials;
	std::map<int, Surface*> adjusted_surfaces;
	std::map<int, Cell*> adjusted_cells;
	std::map<int, Universe*> adjusted_universes;
	std::map<int, Lattice*> adjusted_lattices;

	int uid;


	/**************************************************************************
	 * Ajust the indices for attributes of all cell and lattice
	 * objects in the geometry
	 *************************************************************************/

	/* Adjust the container of surface ids inside each cell to hold the
	 * surfaces' uids */
	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3) {

		Cell* cell = iter3->second;
		int universe = _universes.at(cell->getUniverse())->getUid();

		/* MATERIAL type cells */
		if (cell->getType() == MATERIAL) {
			CellBasic* cell_basic = static_cast<CellBasic*>(cell);
			int material = _materials.at(cell_basic->getMaterial())->getUid();
			cell_basic->adjustKeys(universe, material);
		}

		/* FILL type cells */
		else {
			CellFill* cell_fill = static_cast<CellFill*>(cell);
			int universe_fill = _universes.at(cell_fill->getUniverseFill())->getUid();
			cell_fill->adjustKeys(universe, universe_fill);
		}
	}


	/* Adjust the container of universe ids inside each lattice to hold the
	 * universes' uids */
	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
		Lattice* lattice = iter5->second;
		lattice->adjustKeys();
	}


	/**************************************************************************
	 * Ajust the indices of the containers of geometry objects which are
	 * attributes of this geometry class
	 *************************************************************************/

	/* Adjust material indices to be uids */
	try {
		for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1) {
			uid = iter1->second->getUid();
			Material* material = iter1->second;
			adjusted_materials.insert(std::pair<int, Material*>(uid, material));
		}
		_materials.clear();
		_materials = adjusted_materials;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust material' keys. Backtrace:\n%s", e.what());
	}


	/* Adjust surfaces indices to be uids */
	try {
		for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2) {
			uid = iter2->second->getUid();
			Surface* surface = iter2->second;
			adjusted_surfaces.insert(std::pair<int, Surface*>(uid, surface));
		}
		_surfaces.clear();
		_surfaces = adjusted_surfaces;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust surface' keys. Backtrace:\n%s", e.what());
	}

	/* Adjust cells indices to be uids */
	try {
		for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3) {
			uid = iter3->second->getUid();
			Cell* cell = iter3->second;
			adjusted_cells.insert(std::pair<int, Cell*>(uid, cell));
		}
		_cells.clear();
		_cells = adjusted_cells;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust cell' keys Backtrace:\n%s", e.what());
	}

	/* Adjust universes indices to be uids */
	try {
		for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4) {
			uid = iter4->second->getUid();
			Universe* universe = iter4->second;
			adjusted_universes.insert(std::pair<int, Universe*>(uid, universe));
		}
		_universes.clear();
		_universes = adjusted_universes;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust universes' keys Backtrace:\n%s", e.what());
	}

	/* Adjust lattices indices to be uids */
	try {
		for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
			uid = iter5->second->getUid();
			Lattice* lattice = iter5->second;
			adjusted_universes.insert(std::pair<int, Lattice*>(uid, lattice));
		}
		_lattices.clear();
		_lattices = adjusted_lattices;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust lattices' keys Backtrace:\n%s", e.what());
	}

	return;
}


/**
 * Builds each surfaces list of neighboring cells on the positive and negative side of the
 * surface. This function helps speed up searches for the next cell when for a surface is
 * is crossed while segmenting tracks across the geometry
 */
void Geometry::buildNeighborsLists() {

	log_printf(NORMAL, "Building neighbor cell lists for each surface...\n");

	int count_positive[_surfaces.size()];
	int count_negative[_surfaces.size()];
	std::map<int, Cell*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;

	/* Initialize counts to zero */
	for (int i = 0; i < (int)_surfaces.size(); i++) {
		count_positive[i] = 0;
		count_negative[i] = 0;
	}

	/* Build counts */
	/* Loop over all cells */
	for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {
		std::map<int, Surface*> surfaces = iter1->second->getSurfaces();

		/* Loop over all of this cell's surfaces */
		for (iter2 = surfaces.begin(); iter2 != surfaces.end(); ++iter2) {
			int surface_id = iter2->first;
			bool sense = surface_id > 0;
			surface_id = abs(surface_id);
			if (sense)
				count_positive[surface_id]++;
			else
				count_negative[surface_id]++;
		}
	}

	/* Allocate memory for neighbor lists for each surface */
	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2) {
		int surface = iter2->first;
		if (count_positive[surface] > 0)
			iter2->second->setNeighborPosSize(count_positive[surface]);
		if (count_negative[surface] > 0)
			iter2->second->setNeighborNegSize(count_negative[surface]);
	}

	/* Reinitialize counts to zero */
	for (int i = 0; i < (int)_surfaces.size(); i++) {
		count_positive[i] = 0;
		count_negative[i] = 0;
	}

	/* Loop over all cells */
	for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {
		std::map<int, Surface*> surfaces = iter1->second->getSurfaces();

		/* Loop over all of this cell's surfaces */
		for (iter2 = surfaces.begin(); iter2 != surfaces.end(); ++iter2) {
			int surface = iter2->first;
			bool sense = (surface > 0);
			surface = abs(surface);

			if (sense) {
				count_positive[surface]++;
				iter2->second->setNeighborPos(count_positive[surface], iter1->second);
			}
			else {
				count_negative[surface]++;
				iter2->second->setNeighborNeg(count_negative[surface], iter1->second);
			}
		}
	}

	return;
}


//bool Geometry::findCell(LocalCoords* coords) {
//
//	int universe_id = coords->getUniverse();
//	std::vector<int> cell_ids = _universes.at(universe_id)->getCells();
//	int num_cells = cell_ids.size();
//
//	/* Loop over all cells in this universe */
//	for (int c = 0; c < num_cells; c++) {
//		Cell* cell = _cells.at(cell_ids.at(c));
//
//		if (cellContains(cell, coords)) {
//			/* Set the cell on this level */
//			coords->setCell(cell->getUid());
//
//			/* MATERIAL type cell - lowest level, terminate search for cell */
//			if (cell->getType() == MATERIAL) { }
//
//			/* FILL type cell - cell contains a universe at a lower level
//			 * Update coords to next level and continue search
//			 */
//			else if (cell->getType() == FILL) {
//				LocalCoords* new_coords = new LocalCoords(coords->getX(),coords->getY());
//				coords->setNext(new_coords);
//				coords->setUniverse(cell->getUniverse());
//
//				if (!findCell(coords))
//					return false;
//			}
//
//
//
//			/* UNIVERSE FILL type cell */
//
//
//			/* Lattice ? */
//
//			/* If none of the nested cells contained this coord
//			 * This should not be invoked unless there is a problem with
//			 * the way the geometry is setup
//			 */
//			if (!findCell(coords))
//				return false;
//		}
//	}
//
//	return false;
//
//}


/**
 * Function to determine whether a key already exists in a templated map container
 * @param map the map container
 * @param key the key to check
 * @return true if already in map, false otherwise
 */
template <class K, class V>
bool Geometry::mapContainsKey(std::map<K, V> map, K key) {
	/* Try to access the element at the key */
	try { map.at(key); }

	/* If an exception is thrown, element does not exist */
	catch (std::exception& exc) { return false; }

	/* If no exception is thrown, element does exist */
	return true;
}
