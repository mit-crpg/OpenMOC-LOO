/*
 * Geometry.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Geometry.h"


/**
 * Geometry constructor
 * @param num_sectors the number of angular sectors per cell
 * @param num_rings the number of rings per cell
 * @param sector_offset the angular offset for computing angular sectors
 * @param parser pointer to the parser object
 */
Geometry::Geometry(int num_sectors, int num_rings, double sector_offset, 
		   Parser* parser) {

	_num_sectors = num_sectors;
	_num_rings = num_rings;
	_sector_offset = sector_offset;

	/* Initializing the corners to be infinite  */
	_x_min = 1.0/0.0;
	_y_min = 1.0/0.0;
	_x_max = -1.0/0.0;
	_y_max = -1.0/0.0;

	/* Add each material from parser */
	parser->each_material([this](Material *m) -> void
			 {
				 addMaterial(m);
				 return;
			 });

	/* Add each surface from parser */
	parser->each_surface([this](Surface *s) -> void
			{
				addSurface(s);
				return;
			});//	coords->setLattice(_uid);


	/* Add each cell from parser */
	parser->each_cell([this](Cell *c) -> void
		     {
			     addCell(c);
			     return;
		     });

	/* Add each lattice from parser */
	parser->each_lattice([this](Lattice *l) -> void
			{
				addLattice(l);
				return;
			});

	Universe *univ = _universes.at(0);
	_num_FSRs = univ->computeFSRMaps();
	log_printf(NORMAL, "Number of flat source regions computed: %d", _num_FSRs);
}


/**
 * geometry Destructor clears all memory for materials, surfaces, cells,
 * universes and lattices
 */
Geometry::~Geometry() {
	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;


	for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
		delete iter1->second;
	_materials.clear();

	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
		delete iter2->second;
	_surfaces.clear();

	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
		delete iter3->second;
	_cells.clear();

	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4)
		delete iter4->second;
	_universes.clear();
	_lattices.clear();
}


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
void Geometry::setNumSectors(int num_sectors) {//	coords->setLattice(_uid);

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
 * Returns the number of angular sectors for subdividing flat source regions
 * @param the number of angular sectors
 */
int Geometry::getNumSectors() const {//	coords->setLattice(_uid);

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
 * Returns the number of flat source regions in the geometry
 * @return number of flat source regions
 */
int Geometry::getNumFSRs() const {
	return _num_FSRs;
}


/**
 * Add a material to the geometry
 * @param material a pointer to a material object
 */
void Geometry::addMaterial(Material* material) {

	/* Checks if material with same id has already been added */
	if (_materials.find(material->getId()) != _materials.end())
		log_printf(ERROR, "Cannot add a second material with id = %d",
				material->getId());

	else {
		try {
			_materials.insert(std::pair<int, Material*>(material->getId(),
					material));
			log_printf(INFO, "Added material with id = %d to geometry",
					material->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add material with id = %d. Backtrace:"
					"\n%s", material->getId(), e.what());
		}
	}
}


/**
 * Return a material from the geometry
 * @param id the material id
 * @return a pointer to the material object
 *
 */
Material* Geometry::getMaterial(int id) {
	try {
		return _materials.at(id);
	}
	catch (std::exception & e) {//	coords->setLattice(_uid);

		log_printf(ERROR, "Attempted to retrieve material with id = %d which"
				" does not exist. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a surface to the geometry
 * @param a pointer to the surface object
 */
void Geometry::addSurface(Surface* surface) {

	/* Checks if a surface with the same id has already been added */
	if (_surfaces.find(surface->getId()) != _surfaces.end())
		log_printf(ERROR, "Cannot add a second surface with id = %d",
				surface->getId());

	else {
		try {
			_surfaces.insert(std::pair<int, Surface*>(surface->getId(),
					surface));
			log_printf(INFO, "Added surface with id = %d to geometry",
					surface->getId());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to add surface with id = %d. Backtrace:"
					"\n%s", surface->getId(), e.what());
		}
	}

	/* Use new surface to update the boundaries of the geometry */
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
		log_printf(ERROR, "Attempted to retrieve surface with id = %d which "
				"has not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}


/**
 * Add a cell to the geometry. Checks if the universe the cell is in already
 * exists; if not, it creates one and adds it to the geometry.
 * @param cell a pointer to the cell object
 */
void Geometry::addCell(Cell* cell) {

	/* Checks if a cell with the same id has already been added */
	if (_cells.find(cell->getId()) != _cells.end())
		log_printf(ERROR, "Cannot add a second cell with id = %d",
				cell->getId());

	/* If the cell is filled with a material which does not exist */
	else if (cell->getType() == MATERIAL &&
			_materials.find(static_cast<CellBasic*>(cell)->getMaterial()) ==
						_materials.end()) {

		log_printf(ERROR, "Attempted to add cell with material with id = %d,"
				" but material does not exist",
				static_cast<CellBasic*>(cell)->getMaterial());
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
					"but surface does not exist", iter->first);

		/* The surface does exist, so set the surface pointer in the cell */
		else	int findFSRId(LocalCoords* coords);

			cell->setSurfacePointer(_surfaces.at(surface_id));
	}

	/* Insert the cell into the geometry's cell container */
	try {
		_cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
		log_printf(INFO, "Added cell with id = %d to geometry", cell->getId());
	}
	catch (std::exception &e) {
			log_printf(ERROR, "Unable to add cell with id = %d. Backtrace:"
					"\n%s", cell->getId(), e.what());
	}

	/* Checks if the universe the cell in exists; if not, creates universe */
	if (_universes.find(cell->getUniverse()) == _universes.end()) {
		try {
			Universe* univ = new Universe(cell->getUniverse());
			addUniverse(univ);
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to create a new universe with id = %d "
					"and add it to the geometry. Backtrace:\n%s",
					cell->getUniverse(), e.what());
		}
	}

	/* Adds the cell to the appropriate universe */
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


/*
 * Add a universe to the geometry
 * @param universe a pointer to the universe object
 */
void Geometry::addUniverse(Universe* universe) {
	/* Checks if a universe with the same id has already been added */
	if (_universes.find(universe->getId()) != _universes.end())
		log_printf(ERROR, "Cannot add a second universe with id = %d",
				universe->getId());

	/* Add the universe */
	else {
		try {
			_universes.insert(std::pair<int, Universe*>(universe->getId(),
														universe));
			log_printf(INFO, "Added universe with id = %d to geometry",
					universe->getId());
		}
		catch (std::exception &e) {
				log_printf(ERROR, "Unable to add universe with id = %d. "
						"Backtrace:\n%s", universe->getId(), e.what());
		}
	}
	//	coords->setLattice(_uid);

	/* Checks if any cellfill references this universe and sets its pointer */
	std::map<int, Cell*>::iterator iter;
	for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
		if (iter->second->getType() == FILL) {
			CellFill* cell = static_cast<CellFill*>(iter->second);
			if (cell->getUniverseFillId() == universe->getId())
				cell->setUniverseFillPointer(universe);	int findFSRId(LocalCoords* coords);

		}
	}

	return;
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
 */
void Geometry::addLattice(Lattice* lattice) {
	/* Checks whether a lattice with the same id has already been added */
	if (_lattices.find(lattice->getId()) != _lattices.end())
		log_printf(ERROR, "Cannot add a second lattice with id = %d",
				lattice->getId());

	/* If the universes container already has a universe with the same id */
	else if (_universes.find(lattice->getId()) != _universes.end())
		log_printf(ERROR, "Cannot add a second universe (lattice) with "
				"id = %d", lattice->getId());

	/* Sets the universe pointers for the lattice and c//	coords->setLattice(_uid);
	 * hecks if the lattice
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
			else	int findFSRId(LocalCoords* coords);

				lattice->setUniversePointer(_universes.at(universe_id));
		}
	}

	/* Add the lattice to the geometry's lattices container */
	try {
		_lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
		log_printf(INFO, "Added lattice with id = %d to geometry",
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
 * @param the lattiLocalCoords* ce (universe) id
 * @return a pointer to the lattice object
 */
Lattice* Geometry::getLattice(int id) {
	try {
		return _lattices.at(id);
	}
	catch (std::exception & e) {
		log_printf(ERROR, "Attempted to retrieve lattice with id = %d which "
				"has not been declared. Backtrace:\n%s", id, e.what());
	}
	exit(0);
}
int findFSRId(LocalCoords* coords);


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

	string << "Geometry: width = " << getWidth() << ", height = " <<
			getHeight() << ", base universe id = " << _base_universe <<
			", Bounding Box: (("
			<< _x_min << ", " << _y_min << "), (" << _x_max << ", " << _y_max
			<< ")";

	string << "\n\tMaterials:\n\t\t";
	for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
		string << iter1->second->toString() << "\n\n\t\t";

	string << "\n\tSurfaces:\n\t\t";
	for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
		string << iter2->second->toString() << "\n\t\t";

	string << "\n\tCells:\n\t\t";
	for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
		string << iter3->second->toString() << "\n\t\t";

	string << "\n\tUniverses:\n\t\t";
	for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4) {
		string << iter4->second->toString() << "\n\t\t";
	}

	string << "\n\tLattices:\n\t\t";	int findFSRId(LocalCoords* coords);

	for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
		string << iter5->second->toString()  << "\n\t\t";
	}

	std::string formatted_string = string.str();
	formatted_string.erase(formatted_string.end()-3);

	return formatted_string;
}


/*
 * Prints a string representation of all of the geometry's objects to
 * the console
 */
void Geometry::printString() {
	log_printf(RESULT, "Printing the geometry to the console:\n\t%s",
				toString().c_str());
	return;
}


/* Adjusts the keys for surfaces, cells, universes, and lattices to uids
 */
void Geometry::adjustKeys() {

	log_printf(NORMAL, "Adjusting the keys for the geometry...");

	/* Iterators for all geometry classes */
	std::map<int, Material*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;
	std::map<int, Cell*>::iterator iter3;
	std::map<int, Universe*>::iterator iter4;
	std::map<int, Lattice*>::iterator iter5;

	/* New maps for geometry classes using uids instead of ids as keys */
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
		//	coords->setLattice(_uid);

		/* MATERIAL type cells */
		if (cell->getType() == MATERIAL) {
			CellBasic* cell_basic = static_cast<CellBasic*>(cell);
			int material = _materials.at(cell_basic->getMaterial())->getUid();
			cell_basic->adjustKeys(universe, material);
		}

		/* FILL type cells */
		else {
			CellFill* cell_fill = static_cast<CellFill*>(cell);
			int universe_fill = _universes.at(cell_fill->getUniverseFillId())->
											getUid();
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
	 * Adjust the indices of the containers of geometry objects which are
	 * attributes of this geometry class
	 *************************************************************************/

	/* Adjust material indices to be uids */
	try {
		for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1) {
			uid = iter1->second->getUid();
			Material* material = iter1->second;
			adjusted_materials.insert(std::pair<int, Material*>
														(uid, material));
		}

		/* Reset the materials container to the new map*/
		_materials.clear();
		_materials = adjusted_materials;
	}//	coords->setLattice(_uid);

	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust material' keys. Backtrace:\n%s",
					e.what());
	}


	/* Adjust surfaces indices to be uids */
	try {
		for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2) {
			uid = iter2->second->getUid();
			Surface* surface = iter2->second;
			adjusted_surfaces.insert(std::pair<int, Surface*>(uid, surface));
		}

		/* Reset the surfaces container to the new map*/
		_surfaces.clear();
		_surfaces = adjusted_surfaces;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust surface' keys. Backtrace:\n%s",
					e.what());
	}

	/* Adjust cells indices to be uids */
	try {
		for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3) {
			uid = iter3->second->getUid();
			Cell* cell = iter3->second;
			adjusted_cells.insert(std::pair<int, Cell*>(uid, cell));
		}

		/* Reset the cell container to the new map*/
		_cells.clear();
		_cells = adjusted_cells;	int findFSRId(LocalCoords* coords);

	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust cell' keys Backtrace:\n%s",
					e.what());
	}

	/* Adjust universes indices to be uids */
	try {//	coords->setLattice(_uid);

		for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4) {
			uid = iter4->second->getUid();
			Universe* universe = iter4->second;
			adjusted_universes.insert(std::pair<int, Universe*>
												(uid, universe));
		}

		/* Reset the universes container to the new map*/
		_universes.clear();
		_universes = adjusted_universes;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust universes' keys Backtrace:\n%s",
					e.what());
	}

	/* Adjust lattices indices to be uids */
	try {
		for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5) {
			uid = iter5->second->getUid();
			Lattice* lattice = iter5->second;
			adjusted_lattices.insert(std::pair<int, Lattice*>(uid, lattice));
		}

		/* Reset the lattices container to the new map*/
		_lattices.clear();
		_lattices = adjusted_lattices;
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to adjust lattices' keys Backtrace:\n%s",
					e.what());
	}

	return;
}


/**
 * Builds each surfaces list of neighboring cells on the positive and negative
 * side of the surface. This function helps speed up searches for the next cell
 * when for a surface is is crossed while segmenting tracks across the geometry
 */
void Geometry::buildNeighborsLists() {

	log_printf(NORMAL, "Building neighbor cell lists for each surface...");

	/* Arrays to count the number of surfaces found on the positive/negative
	 * side of each surface */
	int count_positive[_surfaces.size()];
	int count_negative[_surfaces.size()];

	/* Cell and Surface map iterators */
	std::map<int, Cell*>::iterator iter1;
	std::map<int, Surface*>::iterator iter2;

	/* Initialize counts to zero *///	coords->setLattice(_uid);

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
				iter2->second->setNeighborPos(count_positive[surface],
						iter1->second);
			}
			else {
				count_negative[surface]++;
				iter2->second->setNeighborNeg(count_negative[surface],
						iter1->second);
			}
		}
	}

	return;
}


/**
 * Find the cell that this localcoords object is in. This method assumes that
 * the localcoord has coordinates and a universe id. The method will
 * recursively find the localcoord by building a linked list of localcoords
 * from the localcoord passed in as an argument down to the lowest level cell
 * found. In the process it will set the local coordinates for each localcoord
 * in the linked list for the lattice or universe that it is in. If the
 * localcoord is outside the bounds of the geometry or on the boundaries this
 * method will return NULL; otherwise it will return a pointer to the cell
 * that the localcoords is currently in.
 * @param coords pointer to a localcoords object
 * @return returns a pointer to a cell if found, NULL if no cell found
 */
Cell* Geometry::findCell(LocalCoords* coords) {
	int universe_id = coords->getUniverse();
	Universe* univ = _universes.at(universe_id);
	return univ->findCell(coords, _universes);
}


/**
 * Finds the next cell for a localcoords object along a trajectory defined
 * by some angle (in radians from 0 to PI). The method will update the
 * localcoord passed in as an argument to be the one at the boundary of the
 * next cell crossed along the given trajectory. It will do this by
 * recursively building a linked list of localcoords from the localcoord
 * passed in as an argument down to the lowest level cell found. In the
 * process it will set the local coordinates for each localcoord in the
 * linked list for the lattice or universe that it is in. If the
 * localcoord is outside the bounds of the geometry or on the boundaries this
 * method will return NULL; otherwise it will return a pointer to the cell
 * that the localcoords will reach next along its trajectory.
 * @param coords pointer to a localcoords object
 * @param angle the angle of the trajectory
 * @return a pointer to a cell if found, NULL if no cell found
 */
Cell* Geometry::findNextCell(LocalCoords* coords, double angle) {

	Cell* cell = NULL;
	double dist;

	/* Find the current cell */
	cell = findCell(coords);

	/* If the current coords is not in any cell, return NULL */
	if (cell == NULL)
		return NULL;

	/* If the current coords is inside a cell, look for next cell */
	else {
		/* Check the min dist to the next surface in the current cell */
		Point surf_intersection;
		LocalCoords* lowest_level = coords->getLowestLevel();
		dist = cell->minSurfaceDist(lowest_level->getPoint(), angle,
									&surf_intersection);

		/* If the distance returned is not INFINITY, the trajectory will
		 * intersect a surface in the cell */
		if (dist != INFINITY) {
			/* Move LocalCoords just to the next surface in the cell plus an
			 * additional small bit into the next cell */
			coords->updateMostLocal(&surf_intersection);
			double delta_x = cos(angle) * TINY_MOVE;
			double delta_y = sin(angle) * TINY_MOVE;
			coords->adjustCoords(delta_x, delta_y);

			/* Find new cell and return it */
			cell = findCell(coords);
			return cell;
		}

		/* If the distance returned is infinity, the trajectory will not
		 * intersect a surface in the cell. We thus need to readjust to
		 * the localcoord to the base universe and check whether we need
		 * to move to a new lattice cell */
		else if (dist == INFINITY) {

			/* Get the lowest level localcoords in the linked list */
			LocalCoords* curr = coords->getLowestLevel();

			/* Retrace linkedlist from lowest level */
			while (curr != NULL && curr->getUniverse() != 0) {
				curr = curr->getPrev();

				/* If we reach a localcoord in a lattice, delete all lower
				 * level localcoords in linked list and break loop. */
				if (curr->getType() == LAT) {
					curr->prune();
					curr = NULL;
				}
			}

			/* Get the lowest level universe in linkedlist */
			curr = coords->getLowestLevel();

			/* Retrace through the lattices in the localcoord and check for
			 * lattice cell crossings in each one. If we never find a crossing
			 * and reach universe 0 the return NULL since this means we have
			 * reached the edge of the geometry
			 */
			while (curr->getUniverse() != 0) {

				/* If the lowest level localcoords is inside a lattice, find
				 * the next lattice cell */
				if (curr->getType() == LAT) {

					int lattice_id = curr->getLattice();
					Lattice* lattice = _lattices.at(lattice_id);

					cell = lattice->findNextLatticeCell(curr, angle,
															_universes);

					/* If the cell returned is NULL, the localcoords are outside
					 * of the current lattice, so move to a higher level lattice
					 * if there is one */
					if (cell == NULL) {

						/* Delete current lattice */
						curr->getPrev()->prune();

						/* Get the lowest level localcoords in the linked list */
						curr = coords->getLowestLevel();

						/* Retrace linkedlist from lowest level */
						while (curr != NULL && curr->getUniverse() != 0) {
							curr = curr->getPrev();

							/* If we reach a localcoord in a lattice, delete all lower
							 * level localcoords in linked list and break loop. */
							if (curr->getType() == LAT) {
								curr->prune();
								curr = NULL;
							}
						}

						/* Get the lowest level universe in linkedlist */
						curr = coords->getLowestLevel();
					}

					/* If the lowest level universe is not a lattice, then
					 * return the current cell */
					else
						return cell;
				}
			}
		}
	}

	/* If no cell was found, return NULL */
	return NULL;
}


/**
 * This method creates segments within flat source regions in the geometry
 * for a given track. It starts at the beginning of the track and finds
 * successive intersection points with flat source regions as the track passes
 * through the geometry and creates segment structs and adds them to the track
 * @param track a pointer to a track to segmentize
 */
void Geometry::segmentize(Track* track) {


	if (track->getStart()->getX() == 0.0 && track->getStart()->getY() == -2.0 && track->getPhi() == 3*M_PI/4) {

	};

	/* Track starting point coordinates and azimuthal angle */
	double x0 = track->getStart()->getX();
	double y0 = track->getStart()->getY();
	double phi = track->getPhi();

	/* Length of each segment */
	double segment_length;

	/* Use a LocalCoords for the start and end of each segment */
	LocalCoords segment_start(x0, y0);
	LocalCoords segment_end(x0, y0);
	segment_start.setUniverse(0);
	segment_end.setUniverse(0);

	/* Find the cell for the track starting point */
	Cell* curr = findCell(&segment_end);
	Cell* prev;

	/* If starting point was outside the bounds of the geometry */
	if (curr == NULL)
		log_printf(WARNING, "Could not find a cell containing the start point "
				"of this track: %s", track->toString().c_str());

	/* While the segment end localcoords is still within the geometry, move
	 * it to the next cell, create a new segment, and add it to the geometry */
	while (curr != NULL) {

		segment_end.copyCoords(&segment_start);

		/* Find the next cell */
		prev = curr;
		curr = findNextCell(&segment_end, phi);

		/* Find the segment length between the segments start and end points */
		segment_length = segment_end.getPoint()->distance(segment_start.getPoint());

		/* Create a new segment */
		segment* new_segment = new segment;
		new_segment->_length = segment_length;

		new_segment->_region_id = findFSRId(&segment_start);
//		new_segment->_region_id = prev->getUid();
//		new_segment->_region_id = static_cast<CellBasic*>(_cells.at(prev->getId()))->getMaterial();

		/* Checks to make sure that new segment does not have the same start
		 * and end points */
		if (segment_start.getX() == segment_end.getX() &&
				segment_start.getY() == segment_end.getY()) {

			log_printf(ERROR, "Created a segment with the same start and end "
					"point: x = %f, y = %f", segment_start.getX(),
					segment_start.getY());
		}
		/* Update coordinates for start of next segment */
		else {
			log_printf(DEBUG, "Created a new segment with start: x = %f, y = "
					"%f, and end: x = %f, y = %f", segment_start.getX(),
					segment_start.getY(), segment_end.getX(),
					segment_end.getY());
		}

		/* Add the segment to the track */
		track->addSegment(new_segment);
	}

	log_printf(INFO, "Created %d segments for track: %s",
			track->getNumSegments(), track->toString().c_str());

	segment_start.prune();
	segment_end.prune();

	return;
}


int Geometry::findFSRId(LocalCoords* coords) {
	int fsr_id = 0;
	LocalCoords* curr = coords;

	while (curr != NULL) {
		if (curr->getType() == LAT) {
			Lattice* lattice = _lattices.at(curr->getLattice());
			fsr_id += lattice->getFSR(curr->getLatticeX(), curr->getLatticeY());
		}
		else if (curr->getType() == UNIV) {
			Universe* universe = _universes.at(curr->getUniverse());
			fsr_id += universe->getFSR(curr->getCell());
		}
		curr = curr->getNext();
	}

	return fsr_id;
}



/*
 * Function to determine whether a key already exists in a templated map
 * container
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
