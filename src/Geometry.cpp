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
Geometry::Geometry(Parser* parser) {

	_max_seg_length = 0;
	_min_seg_length = INFINITY;

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
			});


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


	/* Generate flat source regions */
	Universe *univ = _universes.at(0);
	_num_FSRs = univ->computeFSRMaps();
	log_printf(INFO, "Number of flat source regions computed: %d", _num_FSRs);

	/* Allocate memory for maps between flat source regions ids and cell or
	 * material ids */
	_FSRs_to_cells = new int[_num_FSRs];
	_FSRs_to_materials = new int[_num_FSRs];

	/* Load maps with cell and material ids */
	for (int r=0; r < _num_FSRs; r++) {
		CellBasic* curr =  static_cast<CellBasic*>
						(findCell(_universes.at(0), r));
		_FSRs_to_cells[r] = curr->getId();
		_FSRs_to_materials[r] = curr->getMaterial();
	}

	_mesh = new Mesh;

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

	delete [] _FSRs_to_cells;
	delete [] _FSRs_to_materials;

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
 * Returns the number of flat source regions in the geometry
 * @return number of flat source regions
 */
int Geometry::getNumFSRs() const {
	return _num_FSRs;
}


/* Return the maximum segment length computed during segmentation
 * return max segment length
 */
double Geometry::getMaxSegmentLength() const {
	return _max_seg_length;
}


/* Return the minimum segment length computed during segmentation
 * return min segment length
 */
double Geometry::getMinSegmentLength() const {
	return _min_seg_length;
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
			/* Check that the sum of the material's absorption and scattering
			 * cross-sections equals its total cross-section */
			//material->checkSigmaT();
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
 * Return an array indexed by FSR ids which contains the corresponding cell ids
 * @return an array map of FSR to cell ids
 */
int* Geometry::getFSRtoCellMap() const {
	return _FSRs_to_cells;
}


/**
 * Return an array indexed by FSR ids which contains the corresponding
 * material ids
 * @return an array map of FSR to material ids
 */
int* Geometry::getFSRtoMaterialMap() const {
	return _FSRs_to_materials;
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
	case VACUUM:
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

	/* Prints error msg if a cell with the same id has already been added */
	if (_cells.find(cell->getId()) != _cells.end())
		log_printf(ERROR, "Cannot add a second cell with id = %d",
				   cell->getId());

	/* Prints error msg if the cell is filled with a non-existing material */
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

		/* Prints error msg if the surface does not exist */
		if (_surfaces.find(surface_id) == _surfaces.end())
			log_printf(ERROR, "Attempted to add cell with surface id = %d, "
					   "but surface does not exist", iter->first);

		/* The surface does exist, so set the surface pointer in the cell */
		else	
			int findFSRId(LocalCoords* coords);

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
			log_printf(INFO, "Created universe = %d", cell->getUniverse());
		}
		catch (std::exception &e) {
			log_printf(ERROR, "Unable to create a new universe with id = %d "
					"and add it to the geometry. Backtrace:\n%s",
					cell->getUniverse(), e.what());
		}
	}

	/* Adds the cell to the appropriate universe */
	_universes.at(cell->getUniverse())->addCell(cell);


/************************* Cutting Up The Cell **********************/
	/* Check if the cell has number of rings; if so, add more cells */
	/* FIXME: need to add error checking */

	static int id = 10000;

	if (cell->getType() == MATERIAL) {

		int t_num_rings = dynamic_cast<CellBasic*>(cell)->getNumRings() + 1;
		if (t_num_rings > 1) {
			log_printf(INFO, "Cell %d has multiple rings; num_rings = %d",
					   cell->getId(), dynamic_cast<CellBasic*>
					   (cell)->getNumRings());
			
			Surface *s;
			CellBasic *c;
			int surface_id, old_id, new_id;
			/* case 1: cell is a circle with one surface */
			if (cell->getNumSurfaces() == 1) {
				double r, r0, r1, rold;
				int i = 2;
				
				/* get cell's radius and compute the radius 
				   of the inner-most circle */
				iter = cells_surfaces.begin();
				surface_id = abs(iter->first);
				r0 =  (dynamic_cast<Circle*>
					   (_surfaces.at(surface_id)))->getRadius();
				r1 = r0 / sqrt(t_num_rings);
				rold = r1;
				
				/* create and add the inner-most circle surface */
				s = new Circle(id, BOUNDARY_NONE, 0.0, 0.0, r1);
				old_id = id;
				addSurface(s);
				log_printf(INFO, "Added new %s", s->toString().c_str());

                /* create and add the inner-most circle cell */
				int tmp = -1 * old_id;
				int *list_surfaces  = &tmp; 
				
				c = new CellBasic
					(old_id, cell->getUniverse(), 1, list_surfaces,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0,
					 dynamic_cast<CellBasic*>(cell)->getNumSectors());
				
				/*
				  c = dynamic_cast<CellBasic*>(cell)->clone
					(old_id, 0, dynamic_cast<CellBasic*>(cell)->getNumSectors());
				c->addSurface(-1*old_id, s);
				*/
				id++;
				addCell(c); /* recursively add more cells if there is any sectors */
				log_printf(INFO, "Added new %s", c->toString().c_str());

				while (i < t_num_rings) {
					/* generate radius for the next circle */
					r = sqrt( rold*rold + ((r0*r0)/t_num_rings) );
	
					/* testing clone */
					//c = dynamic_cast<CellBasic*>(cell)->clone(id);
					//c->addSurface(old_id, s);
				
					/* create and add new surface */
					s = new Circle(id, BOUNDARY_NONE, 0, 0, r);
					addSurface(s);
					new_id = id;

					/* create and add the new ring */				
					int tmp[2] = {old_id, -1 * new_id};
					list_surfaces = &tmp[0];
					c = new CellBasic
						(id, cell->getUniverse(), 2, list_surfaces,
						 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0,
						 dynamic_cast<CellBasic*>(cell)->getNumSectors()); 
					//c->addSurface(-1*new_id, s);
					
					id++;
					addCell(c);
					log_printf(INFO, "Added  %s", c->toString().c_str());	
      
					/* book-keeping */
					rold = r;
					i++;
					old_id = new_id;
				}

				/* update the original circle cell to be the outsidemost cell */
				static_cast<Cell*>(cell)->addSurface(old_id, s); 
				log_printf(INFO, "Update original %s",cell->toString().c_str());
				
			} /* end of case 1*/

			/* case 2: cell is a ring with two surfaces */
			else if (cell->getNumSurfaces() == 2) {
				double r, r01, r02, r1, rold;
				int i = 2, inner_surface, outer_surface;
				
				//std::map<int, Surface*> cells_surfaces = cell->getSurfaces();
	
				/* get cell's two surfaces */
				iter = cells_surfaces.begin();
				int surface_id = abs(iter->first);
				iter++;
				int surface_id2 = abs(iter->first);
			
				/* distinguish which surface is the inner one */
				if (surface_id < surface_id2) {
					inner_surface = surface_id;
					outer_surface = surface_id2;
				}
				else {
					inner_surface = surface_id2;
					outer_surface = surface_id; 
				}
			      
				/* get the cell's two radii */
				r01=(dynamic_cast<Circle*>(_surfaces.at(inner_surface)))
					->getRadius();
				r02=(dynamic_cast<Circle*>(_surfaces.at(outer_surface)))
					->getRadius();
				log_printf(INFO, "Read a ring with radii %f and %f", r01, r02);
			
				/* generate the inner-most radius */
				r1 =  sqrt( r01*r01 + ((r02*r02 - r01*r01)/t_num_rings) );
				rold = r1;

				/* create the inner-most circle surface */
				s = new Circle(id, BOUNDARY_NONE, 0, 0, r1);
				addSurface(s);
				old_id = id;
				log_printf(INFO, "%s", s->toString().c_str());


				/* initialize the inner-most circle cell with the inner radius*/
				int tmp[2] = {-1*old_id, inner_surface};
				int *list_surfaces = &tmp[0];
#if 1
				c = new CellBasic
					(old_id, cell->getUniverse(), 2, list_surfaces, 
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0 ,
					 dynamic_cast<CellBasic*>(cell)->getNumSectors()); 
#else
				c = dynamic_cast<CellBasic*>(cell)->clone
					(old_id, 0, dynamic_cast<CellBasic*>(cell)->getNumSectors());
				c->addSurface(inner_surface, _surfaces.at(inner_surface));
				c->addSurface(-1*old_id, s);
#endif

				id++;
				addCell(c);
				log_printf(INFO, "Added  %s", c->toString().c_str());

				while (i < t_num_rings) {
					/* generate radius for the next circle */
					r = sqrt( rold*rold + ((r02*r02 - r01*r01)/t_num_rings) );

					/* create the new surface and add to the new cell */
					s = new Circle(id, BOUNDARY_NONE, 0, 0, r);
					addSurface(s);	
					new_id = id;

					/* create a new ring cell, and add the old surface before
					   we generate a new one*/
					int tmpp[2] = {old_id, -1 * new_id};
					int *list_s = &tmpp[0];
					c = new CellBasic
						(new_id, cell->getUniverse(), 2, list_s,
						 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0,
						 dynamic_cast<CellBasic*>(cell)->getNumSectors());
					id++;
					addCell(c);
					log_printf(INFO, "Added  %s", c->toString().c_str());	
      
					/* book-keeping */
					rold = r;
					i++;
					old_id = new_id;
				} 

				/* update the original circle cell to be the outside most 
				   ring cell */
				dynamic_cast<Cell*>(cell)->addSurface(old_id, s); 
				log_printf(INFO, "Update original ring %s", cell->toString().c_str());
			} /* end of case 2 */
			/* unsupported surface types */
			else {
				log_printf(ERROR, 
						   "num_rings not supported for these surfaces");	
			}
		} /* end of adding in rings */

		/* begining of adding in sectors */
		int t_num_sectors = dynamic_cast<CellBasic*>(cell)->getNumSectors();
	
		if (t_num_sectors > 0) {
			int *list;
			int num;
			int surface1, surface2, surface3, surface4;
			int surface5, surface6, surface7, surface8;
			Surface *s1, *s2, *s3, *s4, *s5, *s6, *s7, *s8; 
			CellBasic *c1, *c2, *c3, *c4, *c5, *c6, *c7;
			CellBasic *c8, *c9, *c10, *c11, *c12, *c13, *c14, *c15;

			/* generate a list of the current cells */
			std::map<int, Surface*> cells_surfaces = cell->getSurfaces();
			int i = 0;
			num = cell->getNumSurfaces();
			int *tmp = new int[num];
			for (iter = cells_surfaces.begin(); 
				 iter != cells_surfaces.end(); iter++) {
				tmp[i] = iter->first;
				i++;
			}
			list = &tmp[0];

			/* adding in 4 sectors */
			if (t_num_sectors == 4){
				/* generate 2 surfaces */
				surface1 = id;
				log_printf(INFO, "%d", surface1);
				s1 = new Plane(id, BOUNDARY_NONE, 1.0, 1.0, 0);
				addSurface(s1);
				id++;

				surface2 = id;
				log_printf(INFO, "%d", surface2);
				s2 = new Plane(id, BOUNDARY_NONE, 1.0, -1.0, 0);
				addSurface(s2);
				id++;

				/*generate 4 cells */
				c1 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c1);
				c1->addSurface(surface1, s1);
				c1->addSurface(surface2, s2);
				id++;
				
				c2 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c2);
				c2->addSurface(-1*surface1, s1);
				c2->addSurface(-1*surface2, s2);
				id++;
				
				
				c3 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c3);
				c3->addSurface(-1*surface1, s1);
				c3->addSurface(surface2, s2);
				id++;
				
				/* update original cell */
				dynamic_cast<CellBasic*>(cell)->setNumSectors(0);				
				cell->addSurface(surface1, s1);
				cell->addSurface(-1*surface2, s2);
				log_printf(INFO, "original cell is updated to %s",
				cell->toString().c_str());

			} /* end of # sectors = 4 */
			/* adding in 8 sectors */
			else if (t_num_sectors == 8){
				/* generate 4 surfaces */
				surface1 = id;
				log_printf(INFO, "%d", surface1);
				s1 = new Plane(id, BOUNDARY_NONE, 1.0, 1.0, 0);
				addSurface(s1);
				id++;

				surface2 = id;
				log_printf(INFO, "%d", surface2);
				s2 = new Plane(id, BOUNDARY_NONE, 1.0, 0, 0);
				addSurface(s2);
				id++;

				surface3 = id;
				log_printf(INFO, "%d", surface3);
				s3 = new Plane(id, BOUNDARY_NONE, 1.0, -1.0, 0);
				addSurface(s3);
				id++;

				surface4 = id;
				log_printf(INFO, "%d", surface4);
				s4 = new Plane(id, BOUNDARY_NONE, 0.0, 1.0, 0);
				addSurface(s4);
				id++;

				/* generate 7 additional cells */
				c1 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c1);
				c1->addSurface(surface1, s1);
				c1->addSurface(-1 * surface2, s2);
				id++;
				log_printf(INFO, "add cell %s", c1->toString().c_str());
				
				c2 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c2);
				c2->addSurface(surface2, s2);
				c2->addSurface(-1 * surface1, s1);
				id++;
				log_printf(INFO, "add cell %s", c2->toString().c_str());
				
				
				c3 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c3);
				c3->addSurface(surface2, s2);
				c3->addSurface(-1 * surface3, s3);
				id++;
				log_printf(INFO, "add cell %s", c3->toString().c_str());
				
				c4 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c4);
				c4->addSurface(surface3, s3);
				c4->addSurface(-1 * surface2, s2);
				id++;
				log_printf(INFO, "add cell %s", c4->toString().c_str());

				c5 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c5);
				c5->addSurface(surface3, s3);
				c5->addSurface(surface4, s4);
				log_printf(INFO, "add cell %s", c5->toString().c_str());
				id++;

				c6 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c6);
				c6->addSurface(-1 * surface4, s4);
				c6->addSurface(-1 * surface3, s3);
				log_printf(INFO, "add cell %s", c6->toString().c_str());
				id++;

				c7 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c7);
				c7->addSurface(surface4, s4);
				c7->addSurface(-1 * surface1, s1);
				log_printf(INFO, "add cell %s", c7->toString().c_str());
				id++;


				/* update original cell */
				dynamic_cast<CellBasic*>(cell)->setNumSectors(0);				
				cell->addSurface(surface1, s1);
				cell->addSurface(-1 * surface4, s4);
				log_printf(INFO, "original cell is updated to %s",
				cell->toString().c_str()); 
				
			} /* end of # sectors = 8 */
			else if (t_num_sectors == 16){ /* add in # sectors = 16 */
				double pi = 4 * atan(1); 
				double angle = 22.50 / (2 * pi);
				double len = tan(angle);

				/* generate 8 surfaces */
				surface1 = id;
				s1 = new Plane(id, BOUNDARY_NONE, 1.0, 0.0, 0);
				addSurface(s1);
				id++;

				surface2 = id;
				s2 = new Plane(id, BOUNDARY_NONE, 1.0, -len, 0);
				addSurface(s2);
				id++;

				surface3 = id;
				s3 = new Plane(id, BOUNDARY_NONE, 1.0, -1.0, 0);
				addSurface(s3);
				id++;

				surface4 = id;
				s4 = new Plane(id, BOUNDARY_NONE, len, -1.0, 0);
				addSurface(s4);
				id++;

				surface5 = id;
				s5 = new Plane(id, BOUNDARY_NONE, 0, 1.0, 0);
				addSurface(s5);
				id++;

				surface6 = id;
				s6 = new Plane(id, BOUNDARY_NONE, len, 1.0, 0);
				addSurface(s6);
				id++;

				surface7 = id;
				s7 = new Plane(id, BOUNDARY_NONE, 1.0, 1.0, 0);
				addSurface(s7);
				id++;

				surface8 = id;
				s8 = new Plane(id, BOUNDARY_NONE, 1.0, len, 0);
				addSurface(s8);
				id++;

				/* generate 15 additional cells */
				c1 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c1);
				c1->addSurface(surface1, s1);
				c1->addSurface(-1 * surface2, s2);
				id++;
				log_printf(INFO, "add cell %s", c1->toString().c_str());
				
				c2 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c2);
				c2->addSurface(surface2, s2);
				c2->addSurface(-1 * surface1, s1);
				id++;
				log_printf(INFO, "add cell %s", c2->toString().c_str());
				
				
				c3 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c3);
				c3->addSurface(surface2, s2);
				c3->addSurface(-1 * surface3, s3);
				id++;
				log_printf(INFO, "add cell %s", c3->toString().c_str());
				
				c4 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c4);
				c4->addSurface(surface3, s3);
				c4->addSurface(-1 * surface2, s2);
				id++;
				log_printf(INFO, "add cell %s", c4->toString().c_str());

				c5 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c5);
				c5->addSurface(surface3, s3);
				c5->addSurface(-1 * surface4, s4);
				log_printf(INFO, "add cell %s", c5->toString().c_str());
				id++;

				c6 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c6);
				c6->addSurface(surface4, s4);
				c6->addSurface(-1 * surface3, s3);
				log_printf(INFO, "add cell %s", c6->toString().c_str());
				id++;

				c7 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c7);
				c7->addSurface(surface4, s4);
				c7->addSurface(surface5, s5);
				log_printf(INFO, "add cell %s", c7->toString().c_str());
				id++;


				c8 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c8);
				c8->addSurface(-1 * surface5, s5);
				c8->addSurface(-1 * surface4, s4);
				id++;

				c9 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c9);
				c9->addSurface(surface5, s5);
				c9->addSurface(-1 * surface6, s6);
				id++;

				c10 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c10);
				c10->addSurface(surface6, s6);
				c10->addSurface(-1 * surface5, s5);
				id++;


				c11 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c11);
				c11->addSurface(surface6, s6);
				c11->addSurface(-1 * surface7, s7);
				id++;

				c12 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c12);
				c12->addSurface(surface7, s7);
				c12->addSurface(-1 * surface6, s6);
				id++;


				c13 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c13);
				c13->addSurface(surface7, s7);
				c13->addSurface(-1 * surface8, s8);
				id++;

				c14 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c14);
				c14->addSurface(surface8, s8);
				c14->addSurface(-1 * surface7, s7);
				id++;

				c15 = new CellBasic
					(id, cell->getUniverse(), num, list,
					 dynamic_cast<CellBasic*>(cell)->getMaterial(), 0, 0);
				addCell(c15);
				c15->addSurface(surface8, s8);
				c15->addSurface(-1 * surface1, s1);
				id++;

				/* update original cell */
				dynamic_cast<CellBasic*>(cell)->setNumSectors(0);				
				cell->addSurface(surface1, s1);
				cell->addSurface(-1 * surface8, s8);
				log_printf(INFO, "original cell is updated to %s",
				cell->toString().c_str()); 
				
			} /* end of # sectors = 16 */
			/* other number of sectors */
			else {
				log_printf(ERROR,
						   "OpenMOC only supports #sectors = 4, 8, 16.\n"
						   "You entered #sectors = %d", t_num_sectors);
			}

			delete [] tmp;

		} /* end of adding in sections */
		  
	} /* end of material type cell loop */
		
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
						   std::get<0>(lattice->getUniverses().at(i).at(j)));

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

	log_printf(INFO, "Adjusting the keys for the geometry...");

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
	try {

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

	log_printf(INFO, "Building neighbor cell lists for each surface...");

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


/* Find the first cell of a segment with a starting point that is represented
 * by this localcoords object is in. This method assumes that
* the localcoord has coordinates and a universe id. This method will move the
* initial starting point by a small amount along the direction of the track
* in order to ensure that the track starts inside of a distinct FSR rather than
* on the boundary between two of them. The method will recursively find the
* localcoord by building a linked list of localcoords from the localcoord
* passed in as an argument down to the lowest level cell found. In the process
* it will set the local coordinates for each localcoord
* in the linked list for the lattice or universe that it is in.
* @param coords pointer to a localcoords object
* @return returns a pointer to a cell if found, NULL if no cell found
*/
Cell* Geometry::findFirstCell(LocalCoords* coords, double angle) {
	double delta_x = cos(angle) * TINY_MOVE;
	double delta_y = sin(angle) * TINY_MOVE;
	coords->adjustCoords(delta_x, delta_y);
	return findCell(coords);
}


/**
 * Find the cell for an fsr_id. This function calls the recursive function
 * findCell with a pointer to the base level universe 0
 * @param fsr_id a flat source region id
 * @return a pointer to the cell that this fsr is in
 */
Cell* Geometry::findCell(int fsr_id) {
	return findCell(_universes.at(0), fsr_id);
}


/**
 * Find the cell for an fsr_id at a certain universe level. This is a recursive
 * function which is intended to be called with the base universe 0 and an fsr
 * id. It will recursively call itself until it reaches the cell which
 * corresponds to this fsr.
 * @param univ a universe pointer for this fsr's universe level
 * @param fsr_id a flat source region id
 * @return a pointer to the cell that this fsr is in
 */
Cell* Geometry::findCell(Universe* univ, int fsr_id) {

	Cell* cell = NULL;

	/* Check if the FSR id is out of bounds */
	if (fsr_id < -1 || fsr_id > _num_FSRs)
		log_printf(ERROR, "Tried to find the cell for an fsr_id which does not "
							"exist: %d", fsr_id);


	/* If the universe is a SIMPLE type, then find the cell the smallest fsr map
	   entry that is not larger than the fsr_id argument to this function.	*/
	if (univ->getType() == SIMPLE) {
		std::map<int, Cell*>::iterator iter;
		std::map<int, Cell*> cells = univ->getCells();
		Cell* cell_min = NULL;
		int max_id = 0;
		int min_id = INT_MAX;
		int fsr_map_id;

		/* Loop over this universe's cells */
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			fsr_map_id = univ->getFSR(iter->first);
			if (fsr_map_id <= fsr_id && fsr_map_id >= max_id) {
				max_id = fsr_map_id;
				cell = iter->second;
			}
			if (fsr_map_id < min_id) {
				min_id = fsr_map_id;
				cell_min = iter->second;
			}
		}

		/* If the max_id is greater than the fsr_id, there has either been an error
		  or we are at universe 0 and need to go down one level */
		if (max_id > fsr_id) {
			if (cell_min->getType() == MATERIAL)
				log_printf(ERROR, "Could not find cell for fsr_id = %d: "
						"max_id(%d) > fsr_id(%d)", fsr_id, max_id, fsr_id);
			else {
				CellFill* cellfill = static_cast<CellFill*>(cell_min);
				return findCell(cellfill->getUniverseFill(), fsr_id);
			}
		}
		/* Otherwise, decrement the fsr_id and make recursive call to next
		   universe unless an error condition is met */
		else {
			fsr_id -= max_id;
			if (fsr_id == 0 && cell_min->getType() == MATERIAL)
				return cell;
			else if (fsr_id != 0 && cell_min->getType() == MATERIAL)
				log_printf(ERROR, "Could not find cell for fsr_id = %d: "
					"fsr_id = %d and cell type = MATERIAL", fsr_id, fsr_id);
			else {
				CellFill* cellfill = static_cast<CellFill*>(cell_min);
				return findCell(cellfill->getUniverseFill(), fsr_id);
			}
		}
	}

	/* If the universe is a lattice then we find the lattice cell with the
	   smallest fsr map entry that is not larger than the fsr id argument to
	   the function. */
	else {
		Lattice* lat = static_cast<Lattice*>(univ);
		Universe* next_univ = NULL;
		int num_y = lat->getNumY();
		int num_x = lat->getNumX();
		int max_id = 0;
		int fsr_map_id;

		/* Loop over all lattice cells */
		for (int i = 0; i < num_y; i++) {
			for (int j = 0; j < num_x; j++) {
				fsr_map_id = lat->getFSR(j, i);
				if (fsr_map_id <= fsr_id && fsr_map_id >= max_id) {
					max_id = fsr_map_id;
					next_univ = lat->getUniverse(j, i);
				}
			}
		}

		/* If the max_id is out of bounds, then query failed */
		if (max_id > fsr_id || next_univ == NULL)
			log_printf(ERROR, "No lattice cell found for fsr = %d, max_id = "
						"%d", fsr_id, max_id);

		/* Otherwise update fsr_id and make recursive call to next level */
		fsr_id -= max_id;
		return findCell(next_univ, fsr_id);
	}

	return cell;
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
			LocalCoords test(0,0);

			/* Move LocalCoords just to the next surface in the cell plus an
			 * additional small bit into the next cell */
			double delta_x = cos(angle) * TINY_MOVE;
			double delta_y = sin(angle) * TINY_MOVE;

			/* Copy coords to the test coords before moving it by delta and
			 * finding the new cell it is in - do this for testing purposes
			 * in case the new cell found is NULL or is in a new lattice cell*/
			coords->copyCoords(&test);
			coords->updateMostLocal(&surf_intersection);
			coords->adjustCoords(delta_x, delta_y);

			/* Find new cell and return it */
			cell = findCell(coords);

			/* Check if cell is null - this means that intersection point
			 * is outside the bounds of the geometry and the old coords
			 * should be restored so that we can look for the next
			 * lattice cell */
			LocalCoords* test_curr = test.getLowestLevel();
			LocalCoords* coords_curr = coords->getLowestLevel();

			while (test_curr != NULL && test_curr->getUniverse() != 0 &&
					coords_curr != NULL && coords_curr->getUniverse() !=0){
					/* Check if the next cell found is in the same lattice cell
					 * as the previous cell */
				if (coords_curr->getType() == LAT &&
						test_curr->getType() == LAT) {
					if (coords_curr->getLatticeX() != test_curr->getLatticeX()
					|| coords_curr->getLatticeY() != test_curr->getLatticeY())
						dist = INFINITY;
						break;
				}

				test_curr = test_curr->getPrev();
				coords_curr = coords_curr->getPrev();
			}

			/* If the cell is null then we should reset and find the next lattice
			 * cell rather than return this cell */
			if (cell == NULL)
				dist = INFINITY;

			/* If the distance is not INFINITY then the new cell found is the one
			 * to return
			 */
			if (dist != INFINITY) {
				test.prune();
				return cell;
			}

			/* If the distance is not INFINITY then the new cell found is not
			 * the one to return and we should move to a new lattice cell */
			else
				test.copyCoords(coords);

			test.prune();
		}

		/* If the distance returned is infinity, the trajectory will not
		 * intersect a surface in the cell. We thus need to readjust to
		 * the localcoord to the base universe and check whether we need
		 * to move to a new lattice cell */
		if (dist == INFINITY) {

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
	Cell* curr = findFirstCell(&segment_end, phi);
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
		segment_length = segment_end.getPoint()
			->distance(segment_start.getPoint());

		/* Create a new segment */
		segment* new_segment = new segment;
		new_segment->_length = segment_length;
		new_segment->_material = _materials.at(static_cast<CellBasic*>(prev)
											   ->getMaterial());

		/* Update the max and min segment lengths */
		if (segment_length > _max_seg_length)
			_max_seg_length = segment_length;
		if (segment_length < _min_seg_length)
			_min_seg_length = segment_length;

		log_printf(DEBUG, "segment start x = %f, y = %f, end x = %f, y = %f",
				segment_start.getX(), segment_start.getY(), segment_end.getX(),
				segment_end.getY());

		new_segment->_region_id = findFSRId(&segment_start);

		/* Gets pointer to mesh surfaces that the segment crosses */
		if (_run_cmfd || _run_loo){
			new_segment->_mesh_surface_fwd = 
				_mesh->findMeshSurface(new_segment->_region_id, &segment_end);
			new_segment->_mesh_surface_bwd = 
				_mesh->findMeshSurface(new_segment->_region_id, &segment_start);
		}

		/* Checks to make sure that new segment does not have the same start
		 * and end points */
		if (segment_start.getX() == segment_end.getX() &&
				segment_start.getY() == segment_end.getY()) {

			log_printf(ERROR, "Created a segment with the same start and end "
					"point: x = %f, y = %f", segment_start.getX(),
					segment_start.getY());
		}

		/* Add the segment to the track */
		track->addSegment(new_segment);
	}

	log_printf(INFO, "Created %d segments for track: %s",
			track->getNumSegments(), track->toString().c_str());

	segment_start.prune();
	segment_end.prune();

	log_printf(DEBUG, "max segment length: %f", _max_seg_length);
	log_printf(DEBUG, "min segment length: %f", _min_seg_length);

	return;
}


/**
 * Find and return the id of the flat source region that this localcoords
 * object is inside of
 * @param coords a localcoords object returned from the findCell method
 */
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


/**
 * This function calls the compress cross-sections method for each
 * Material in the Geometry. These methods will find the minimum
 * and maximum indices for the non-zero elements of each cross-section
 * so that fixed source iteration does not waste time on these elements.
 */
void Geometry::compressCrossSections() {

	log_printf(NORMAL, "Compressing cross-sections...");

	std::map<int, Material*>::iterator iter;

	/* Loop over materials */
	for (iter = _materials.begin(); iter != _materials.end(); ++iter)
		iter->second->compressCrossSections();

	return;
}


/**
 * This method is called from the Solver after fixed source iteration
 * to compute the powers (fission rates) for each lattice cell (ie, the pin
 * and assembly powers for most geometries). The method stores the pin powers
 * mapped by FSR id in the second parameter, FSRs_to_pin_powers
 * @param FSRs_to_powers an array of the fission rate inside a given FSR
 * @param FSRs_to_pin_powers an array of the fission rate of the lattice cell
 * this FSR is within
 */
void Geometry::computePinPowers(double* FSRs_to_powers,
								double* FSRs_to_pin_powers) {

	/* Get the base universe */
	Universe* univ = _universes.at(0);

	/* Create a file prefix for the output files to store all the pin powers */
	std::string file_prefix = "PinPowers/universe0";

	/* Make call to recursive function to compute powers at each
	 * level of lattice */
	computePinPowers(univ, (char*)file_prefix.c_str(), 0, FSRs_to_powers,
														FSRs_to_pin_powers);

	return;
}


/**
 * This is a recursive function which computes the powers of all of the FSRs
 * inside a given universe. This function handles both lattices and regular
 * type universes and saves the powers computed for each lattice cell in a
 * file.
 * @param univ a pointer to the universe of interest
 * @param output_file_prefix the prefix for the output file to save the powers
 * @param FSR_id the FSR id prefix from the previous level's FSR map
 * @param FSRs_to_powers array of the fission rates for each FSR
 * @param FSRs_to_pin_powers array of the fission rates for the lattice cell
 * that each FSR is within
 */
double Geometry::computePinPowers(Universe* univ, char* output_file_prefix,
		int FSR_id, double* FSRs_to_powers, double* FSRs_to_pin_powers) {

	/* Power starts at 0 and is incremented for each FSR in this universe */
	double power = 0;

	bool non_zero_power;

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE) {
		std::map<int, Cell*> cells = univ->getCells();
		std::map<int, int> _region_map;
		std::vector<int> fsr_ids;
		Cell* curr;

		/* For each of the cells inside the lattice, check if it is
		 * material or fill type */
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* If the current cell is a MATERIAL type cell, pull its
			 * FSR id from the fsr map and increment the power by the
			 * power for that FSR
			 */
			if (curr->getType() == MATERIAL) {
				int fsr_id = univ->getFSR(curr->getId()) + FSR_id;
				fsr_ids.push_back(fsr_id);
				power += FSRs_to_powers[fsr_id];
			}

			/* If the current cell is a FILL type cell, pull its
			 * FSR id from the fsr map
			 */
			else {
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				int fsr_id = univ->getFSR(curr->getId()) + FSR_id;

				power += computePinPowers(universe_fill, output_file_prefix, fsr_id,
										FSRs_to_powers, FSRs_to_pin_powers);
			}
		}

		/* Loop over all of the FSR ids stored for MATERIAL type cells
		 * and save their pin powers in the FSRs_to_pin_powers map */
		for (int i=0; i < (int)fsr_ids.size(); i++) {
			int fsr_id = fsr_ids.at(i);
			FSRs_to_pin_powers[fsr_id] = power;
		}
	}

	/* If the universe is a LATTICE type universe */
	else {
		Lattice* lattice = static_cast<Lattice*>(univ);
		Universe* curr;
		int num_x = lattice->getNumX();
		int num_y = lattice->getNumY();
		int fsr_id;
		double cell_power = 0;

		/* Create an output file to write this lattice's pin powers to within
		 * a new directory called PinPowers */
		mkdir("PinPowers", S_IRWXU);
		std::stringstream output_file_name;
		output_file_name << output_file_prefix <<
				"_lattice" << lattice->getId() << "_power.txt";
		FILE* output_file = fopen(output_file_name.str().c_str(), "w");

		non_zero_power = false;

		/* Loop over all lattice cells in this lattice */
		for (int i = num_y-1; i > -1; i--) {
			for (int j = 0; j < num_x; j++) {

				/* Get a pointer to the current lattice cell */
				curr = lattice->getUniverse(j, i);

				/* Get the FSR id prefix for this lattice cell */
				fsr_id = lattice->getFSR(j, i) + FSR_id;

				/* Create an output filename for this cell's power */
				std::stringstream file_prefix;
				file_prefix << output_file_prefix << "_lattice" <<
						lattice->getId() << "_x" << j << "_y" << i;

				/* Find this lattice cell's power */
				cell_power = computePinPowers(curr,
								(char*)file_prefix.str().c_str(),
								fsr_id, FSRs_to_powers, FSRs_to_pin_powers);

				/* Write this lattice cell's power to the output file */
				fprintf(output_file, "%f, ", cell_power);

				power += cell_power;

				/* Check if a nonzero power has been computed */
				if (power > 0.0)
					non_zero_power = true;
			}
			/* Move to the next line in the output file */
			fprintf(output_file, "\n");
		}

		fclose(output_file);

		/* Delete this output file if none of the powers were nonzero */
		if (!non_zero_power)
			remove(output_file_name.str().c_str());

	}

	return power;
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


/**
 * This is a recursive function which makes a mesh for solving the
 * Course Mesh Finite Difference (CMFD) diffusion equations. The CMFD mesh must
 * be a structured mesh and defined as a certain lattice level (counting from
 * the top). This function takes in a pointer to a Mesh object and defines
 * the cell dimensions (cellWidth and CellHeight) and geometric dimensions
 * (width and height) of the mesh. This function adds new MeshCell objects
 * to the Mesh and defines the values in each MeshCell.
 * @param mesh a pointer to the Mesh object
 */
void Geometry::makeCMFDMesh(Mesh* mesh, int numAzim, bool multigroup, bool printMatrices, int cmfdLevel){
	log_printf(NORMAL, "Making mesh at level %i...", cmfdLevel);

	int max_cmfd_level = 0;
	Universe* univ = _universes.at(0);

	/* find the mesh depth of the geometry */
	max_cmfd_level = findMeshDepth(univ, max_cmfd_level);
	log_printf(INFO, "Maximum cmfd mesh depth is: %i level(s)", max_cmfd_level);

	/* if cmfd_level > CMFD_LEVEL throw errror message */
	if (cmfdLevel > max_cmfd_level){
		log_printf(ERROR, "Specified CMFD LEVEL (%i) is greater than the "
				"maximum mesh depth (%i). Please reset CMFD_LEVEL to be <= %i.",
				cmfdLevel, max_cmfd_level, max_cmfd_level);
	}

	/* find cell width and height at CMFD_LEVEL lattice */
	int width = 0;
	int height = 0;
	findMeshHeight(univ, &height, cmfdLevel);
	univ = _universes.at(0);
	findMeshWidth(univ, &width, cmfdLevel);

	/* set mesh boundary conditions */
	mesh->setBoundary(getSurface(1)->getBoundary(), 0);
	mesh->setBoundary(getSurface(2)->getBoundary(), 2);
	mesh->setBoundary(getSurface(3)->getBoundary(), 1);
	mesh->setBoundary(getSurface(4)->getBoundary(), 3);

	/* set the cell and geometric width and height of mesh */
	mesh->setCellHeight(height);
	mesh->setCellWidth(width);
	mesh->setHeight(getHeight());
	mesh->setWidth(getWidth());
	mesh->makeMeshCells();

	log_printf(NORMAL, "Made mesh with %i mesh cells...", height*width);

	log_printf(DEBUG, "mesh cell width: %i", _mesh->getCellWidth());
	log_printf(DEBUG, "mesh cell height: %i", _mesh->getCellHeight());

	/* make a vector of FSR ids in each mesh cell */
	int meshCellNum = 0;
	log_printf(DEBUG, "defining mesh...");
	defineMesh(mesh, univ, cmfdLevel, &meshCellNum, 0, true, 0);
	mesh->setCellBounds();
	mesh->setFSRBounds(getSurface(1)->getBoundary(), 
					   getSurface(2)->getBoundary(), 
					   getSurface(3)->getBoundary(), 
					   getSurface(4)->getBoundary());

	log_printf(DEBUG, "Setting multigroup and print flags");
	mesh->setMultigroup(multigroup);
	mesh->setPrintMatrices(printMatrices);

	return;
}


/**
 * This is a recursive function which stores the fsr_id's of all flat source 
 * regions located in a MeshCell object in a vector. The vector is owned by 
 * the MeshCell object.
 * @param univ a pointer to a the universe that contains the MeshCell
 * @param meshCell a pointer to the MeshCell object
 * @param fsr_id a pointer to an integer that is set to first fsr_id in this
 * universe
 */
void Geometry::findFSRs(Universe* univ, MeshCell* meshCell, int *fsr_id){

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE) {
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;

		/* For each of the cells inside the lattice, check if it is
		 * material or fill type */
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* If the current cell is a MATERIAL type cell, store its fsr_id */
			if (curr->getType() == MATERIAL) {
				log_printf(DEBUG, "pushing back fsr id: %i", *fsr_id);
				meshCell->addFSR(*fsr_id);
				log_printf(DEBUG, "size of fsr list: %i", 
						   (int) meshCell->getFSRs()->size());
				*fsr_id += 1;
			}

			/* If the current cell is a FILL type cell recursively call 
			 * findFSRs */
			else {
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				findFSRs(universe_fill, meshCell, fsr_id);
			}
		}
	}

	/* If the universe is a LATTICE type universe recursively call findFSRs */
	else {
		Lattice* lattice = static_cast<Lattice*>(univ);
		Universe* curr;
		int num_x = lattice->getNumX();
		int num_y = lattice->getNumY();
		int baseFSR = *fsr_id;

		/* Loop over all lattice cells in this lattice */
		for (int i = num_y-1; i > -1; i--) {
			for (int j = 0; j < num_x; j++) {

				/* Get a pointer to the current lattice cell */
				curr = lattice->getUniverse(j, i);
				log_printf(DEBUG, "getting lattice fsr: %i", 
						   lattice->getFSR(j,i));
				*fsr_id = baseFSR + lattice->getFSR(j,i);

				/* find all FSRs in this lattice */
				findFSRs(curr, meshCell, fsr_id);
			}
		}
	}
}

/**
 * This is a recursive function which defines all the parameters of the MeshCell objects
 * in a Mesh. This function takes in the base universe (universe 0) and descends until
 * it reaches the universe level of the CMFD mesh. Then, the function loops over all
 * the cells in Lattice and defines the corresponding MeshCell object for each Lattice
 * Cell.
 * @param univ a pointer to a the base universe (universe 0)
 * @param depth the number of lattices that must be descended to reach the CMFD mesh level.
 * @param meshCellNum a pointer to an integer used to store the index of the current
 * MeshCell object.
 * @param row the current row of the parent lattice
 * @param base bool indicating whether the current lattice is the highest level lattice
 * @param fsr_id a pointer to an integer that is set to first fsr_id in this
 * universe
 */
/* define the MeshCell objects (_width, _height, and _FRSs)  */
void Geometry::defineMesh(Mesh* mesh, Universe* univ, int depth, int* meshCellNum, int row, bool base, int fsr_id){

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE){
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;

		/* For each of the cells inside the lattice, check if it is
		 * material or fill type */
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;CellFill* fill_cell = static_cast<CellFill*>(curr);
			Universe* universe_fill = fill_cell->getUniverseFill();
			defineMesh(mesh, universe_fill, depth, meshCellNum, row, base, fsr_id);
		}
	}

	/* If the universe is a LATTICE type universe */
	else {

		Lattice* lattice = static_cast<Lattice*>(univ);
		Universe* curr;
		int num_x = lattice->getNumX();
		int num_y = lattice->getNumY();
		log_printf(DEBUG, "numx: %i numy: %i", num_x, num_y);

		/* if the current LATTICE is the CMFD mesh lattice */
		if (depth == 1){
			/* if the current LATTICE is the base lattice */
			if (base == true){
				for (int i = num_y-1; i > -1; i--) {
					for (int j = 0; j < num_x; j++) {
						curr = lattice->getUniverse(j, i);
						fsr_id = lattice->getFSR(j,i);
						log_printf(DEBUG, "added FSR id to counter -> fsr id:"
								   " %i", fsr_id);

						/* store the fsr_ids of the FSRs in this LATTICE 
						 * in a MeshCell object */
						findFSRs(curr, mesh->getCells(*meshCellNum), &fsr_id);
						double w = lattice->getWidthX();
						double h = lattice->getWidthY();
						double l = sqrt(w * w + h * h) / 2;
						l /= 0.798184;
						mesh->getCells(*meshCellNum)->setWidth(w);
						mesh->getCells(*meshCellNum)->setHeight(h);
						mesh->getCells(*meshCellNum)->setL(l);
						log_printf(DEBUG, "mesh cell: %i, width: %f,"
								   " height: %f", *meshCellNum, w, h);
						
						/* increments the # of mesh cells */
						*meshCellNum = *meshCellNum + 1;
					}
				}
			}
			/* if the current LATTICE is not the base lattice */
			else{
				int baseFSR = fsr_id;
				for (int j = 0; j < num_x; j++) {
					curr = lattice->getUniverse(j, row);
					fsr_id = baseFSR + lattice->getFSR(j,row);
					log_printf(DEBUG, "set fsr id to: %i", fsr_id);

					/* store the fsr_ids of the FSRs in this LATTICE 
					 * in a MeshCell object */
					findFSRs(curr, mesh->getCells(*meshCellNum), &fsr_id);
					double w = lattice->getWidthX();
					double h = lattice->getWidthY();
					double l = 0.5 * sqrt(w * w + h * h) / SIN_THETA_45;	
					mesh->getCells(*meshCellNum)->setWidth(w);
					mesh->getCells(*meshCellNum)->setHeight(h);
					mesh->getCells(*meshCellNum)->setL(l);

					log_printf(DEBUG, "mesh cell: %i, width: %f,"
							   " height: %f", *meshCellNum, w, h);
					
					/* increments the # of mesh cells */
					*meshCellNum = *meshCellNum + 1;
				}
			}
		}

		/* if the current LATTICE is not the CMFD mesh lattice */
		else {
			base = false;
			for (int i = num_y-1; i > -1; i--) {
				curr = lattice->getUniverse(0, i);
				int nextHeight = nextLatticeHeight(curr);
				log_printf(DEBUG, "next height: %i", nextHeight);
				for (int k = nextHeight-1; k > -1; k--) {
					for (int j = 0; j < num_x; j++) {
						curr = lattice->getUniverse(j, i);
						fsr_id = lattice->getFSR(j,i);

						/* recursively call defineMesh until LATTICE level of CMFD mesh is reached */
						defineMesh(mesh, curr, depth - 1, meshCellNum, k, base, fsr_id);
					}
				}
			}
		}
	}

	return;
}


/**
 * This is a recursive function that finds the cellHeight of the next lowest LATTICE in the
 * given universe.
 * @param univ a pointer to a the current universe
 */
int Geometry::nextLatticeHeight(Universe* univ){

	int height = 1;

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE){
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;
		std::map<int, Cell*>::iterator iter;
		iter = cells.begin();
		curr = iter->second;

		/* IF the cell is FILL type recursively call nextLatticeHeight */
		if (curr->getType() == FILL){
			CellFill* fill_cell = static_cast<CellFill*>(curr);
			Universe* universe_fill = fill_cell->getUniverseFill();
			height = nextLatticeHeight(universe_fill);
		}
	}

	/* If the universe is a LATTICE type universe return the height */
	else {
		Lattice* lattice = static_cast<Lattice*>(univ);
		height = lattice->getNumY();
	}

	return height;
}



/**
 * This is a recursive function that finds the cellHeight of the LATTICE at
 * the CMFD mesh level
 * @param univ a pointer to a the base universe (universe 0)
 * @param height a pointer to the accumulator for the Mesh height
 * @param depth the number of lattices that must be descended to reach the CMFD mesh level.
 */
void Geometry::findMeshHeight(Universe* univ, int* height, int depth){

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE){
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* IF the cell is FILL type recursively call findMeshHeight */
			if (curr->getType() == FILL){
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				findMeshHeight(universe_fill, height, depth);
			}
		}
	}

	/* If the universe is a LATTICE type universe */
	else {
		/* check to see if CMFD mesh lattice depth has been reached */
		if (depth > 0){

			Lattice* lattice = static_cast<Lattice*>(univ);
			Universe* curr;
			int num_y = lattice->getNumY();

			/* if the current LATTICE is the CMFD mesh lattice add its numY to height accumulator */
			if (depth == 1){
				*height = *height + num_y;
			}
			else{
				depth = depth - 1;
				for (int i = num_y-1; i > -1; i--) {
					curr = lattice->getUniverse(0, i);

					/* find the height of the current universe */
					findMeshHeight(curr, height, depth);
				}
			}
		}
	}

	return;
}


/**
 * This is a recursive function that finds the cellWidth of the LATTICE at
 * the CMFD mesh level.
 * @param univ a pointer to a the base universe (universe 0)
 * @param width a pointer to the accumulator for the Mesh width
 * @param depth the number of lattices that must be descended to reach the CMFD mesh level.
 */
void Geometry::findMeshWidth(Universe* univ, int* width, int depth){

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE){
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* IF the cell is FILL type recursively call findMeshWidth */
			if (curr->getType() == FILL){
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				findMeshWidth(universe_fill, width, depth);
			}
		}
	}

	/* If the universe is a LATTICE type universe */
	else {
		/* check to see if CMFD mesh lattice depth has been reached */
		if (depth > 0){

			Lattice* lattice = static_cast<Lattice*>(univ);
			Universe* curr;
			int num_x = lattice->getNumX();
			int num_y = lattice->getNumY();
			int i = num_y-1;

			/* if the current LATTICE is the CMFD mesh lattice add its numX to width accumulator */
			if (depth == 1){
				*width = *width + num_x;
			}
			else{
				depth = depth - 1;
				for (int j = 0; j < num_x; j++) {
					curr = lattice->getUniverse(j, i);

					/* find the width of the current universe */
					findMeshWidth(curr, width, depth);
				}
			}
		}
	}

	return;
}

/**
 * This is a recursive function that finds the depth of the geometry mesh.
 * @param univ a pointer to a the base universe (universe 0)
 * @param cmfd_level a pointer to the accumulator for the Mesh level
 */
int Geometry::findMeshDepth(Universe* univ, int cmfd_level){

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE){
		std::map<int, Cell*> cells = univ->getCells();
		Cell* curr;
		std::map<int, Cell*>::iterator iter;

		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* IF the cell is FILL type recursively call findMeshWidth */
			if (curr->getType() == FILL){
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				cmfd_level = findMeshDepth(universe_fill, cmfd_level);
			}
		}
	}

	/* If the universe is a LATTICE type universe */
	else {
		Lattice* lattice = static_cast<Lattice*>(univ);
		Universe* curr;
		int num_x = lattice->getNumX();
		int num_y = lattice->getNumY();
		cmfd_level++;
		int i = num_y-1;
		int j = num_x-1;
		int min = 0;
		int levels = 0;

		/* loop through lattice cells and recursively call  */
		for (int x = 0; x <= j; x++){
			for (int y = 0; y <= i; y++){

				curr = lattice->getUniverse(x, y);

				/* find the width of the current universe */
				levels = findMeshDepth(curr, cmfd_level);

				if (x == 0 && y == 0){
					min = levels;
				}
				else{
					min = std::min(min, levels);
				}

			}
		}

		cmfd_level = min;
	}

	return cmfd_level;
}


Mesh* Geometry::getMesh(){
	return _mesh;
}


/**
 * This method is called from the Solver after fixed source iteration
 * to compute the powers (fission rates) for each lattice cell (ie, the pin
 * and assembly powers for most geometries). The method stores the pin powers
 * mapped by FSR id in the second parameter, FSRs_to_pin_powers
 * @param FSRs_to_powers an array of the fission rate inside a given FSR
 * @param FSRs_to_pin_powers an array of the fission rate of the lattice cell
 * this FSR is within
 */
void Geometry::computePinAbsorption
(double* FSRs_to_absorption[NUM_ENERGY_GROUPS + 1],
 double* FSRs_to_pin_absorption[NUM_ENERGY_GROUPS + 1]) {

	/* Get the base universe */
	Universe* univ = _universes.at(0);

	/* Create a file prefix for the output files to store all the pin powers */
	std::string file_prefix = "SigmaA/universe0";

	/* Make call to recursive function to compute powers at each
	 * level of lattice */
	computePinAbsorption(univ, (char*)file_prefix.c_str(), 0, 
						 FSRs_to_absorption, FSRs_to_pin_absorption);

	return;
}


/**
 * This is a recursive function which computes the powers of all of the FSRs
 * inside a given universe. This function handles both lattices and regular
 * type universes and saves the powers computed for each lattice cell in a
 * file.
 * @param univ a pointer to the universe of interest
 * @param output_file_prefix the prefix for the output file to save the powers
 * @param FSR_id the FSR id prefix from the previous level's FSR map
 * @param FSRs_to_powers array of the fission rates for each FSR
 * @param FSRs_to_pin_powers array of the fission rates for the lattice cell
 * that each FSR is within
 */
double Geometry::computePinAbsorption
(Universe* univ, char* output_file_prefix, int FSR_id, 
 double* FSRs_to_absorption[NUM_ENERGY_GROUPS + 1], 
 double* FSRs_to_pin_absorption[NUM_ENERGY_GROUPS + 1]) {

	/* Power starts at 0 and is incremented for each FSR in this universe */
	double sigma_a = 0;

	/* If the universe is a SIMPLE type universe */
	if (univ->getType() == SIMPLE) {
		std::map<int, Cell*> cells = univ->getCells();
		std::map<int, int> _region_map;
		std::vector<int> fsr_ids;
		Cell* curr;

		/* For each of the cells inside the lattice, check if it is
		 * material or fill type */
		std::map<int, Cell*>::iterator iter;
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			curr = iter->second;

			/* If the current cell is a MATERIAL type cell, pull its
			 * FSR id from the fsr map and increment the power by the
			 * power for that FSR
			 */
			if (curr->getType() == MATERIAL) {
				int fsr_id = univ->getFSR(curr->getId()) + FSR_id;
				fsr_ids.push_back(fsr_id);
				sigma_a += FSRs_to_absorption[0][fsr_id];
			}

			/* If the current cell is a FILL type cell, pull its
			 * FSR id from the fsr map
			 */
			else {
				CellFill* fill_cell = static_cast<CellFill*>(curr);
				Universe* universe_fill = fill_cell->getUniverseFill();
				int fsr_id = univ->getFSR(curr->getId()) + FSR_id;

				sigma_a += computePinAbsorption
					(universe_fill, output_file_prefix, fsr_id,
					 FSRs_to_absorption, FSRs_to_pin_absorption);
			}
		}

		/* Loop over all of the FSR ids stored for MATERIAL type cells
		 * and save their pin powers in the FSRs_to_pin_powers map */
		for (int i=0; i < (int)fsr_ids.size(); i++) {
			int fsr_id = fsr_ids.at(i);
			FSRs_to_pin_absorption[0][fsr_id] = sigma_a;
		}
	}

	/* If the universe is a LATTICE type universe */
	else {
		Lattice* lattice = static_cast<Lattice*>(univ);
		Universe* curr;
		int num_x = lattice->getNumX();
		int num_y = lattice->getNumY();
		int fsr_id;
		double cell_absorption = 0;

		/* Create an output file to write this lattice's pin powers to within
		 * a new directory called PinPowers */
		mkdir("SigmaA", S_IRWXU);
		std::stringstream output_file_name;
		output_file_name << output_file_prefix <<
				"_lattice" << lattice->getId() << "_sigmaA.txt";
		FILE* output_file = fopen(output_file_name.str().c_str(), "w");

		/* Loop over all lattice cells in this lattice */
		for (int i = num_y-1; i > -1; i--) {
			for (int j = 0; j < num_x; j++) {

				/* Get a pointer to the current lattice cell */
				curr = lattice->getUniverse(j, i);

				/* Get the FSR id prefix for this lattice cell */
				fsr_id = lattice->getFSR(j, i) + FSR_id;

				/* Create an output filename for this cell's power */
				std::stringstream file_prefix;
				file_prefix << output_file_prefix << "_lattice" <<
						lattice->getId() << "_x" << j << "_y" << i;

				/* Find this lattice cell's power */
				cell_absorption = computePinAbsorption
					(curr, (char*)file_prefix.str().c_str(),
					 fsr_id, FSRs_to_absorption, FSRs_to_pin_absorption);

				/* Write this lattice cell's power to the output file */
				fprintf(output_file, "%f, ", cell_absorption);

				sigma_a += cell_absorption;
			}
			/* Move to the next line in the output file */
			fprintf(output_file, "\n");
		}

		fclose(output_file);
	}

	return sigma_a;
}

void Geometry::setCmfd(bool runCmfd){
	_run_cmfd = runCmfd;
}

void Geometry::setLoo(bool runLoo){
	_run_loo = runLoo;
}


