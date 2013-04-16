/*
 * Cmfd.cpp
 *
 *  Created on: September 13, 2012
 *      Author: Sam Shaner
 *				MIT, Course 22
 *              shaner@mit.edu
 */

#include "Cmfd.h"

/**
 * Transient constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the trackgenerator
 */
Cmfd::Cmfd(Geometry* geom, Plotter* plotter, Mesh* mesh, bool updateFlux, bool runCmfd) {

	_geom = geom;
	_plotter = plotter;
	_mesh = mesh;
	_update_flux = updateFlux;
	_quad = new Quadrature(TABUCHI);

	if (runCmfd){

		PetscInt size1, size2;
		int ng = NUM_ENERGY_GROUPS;
		if (_mesh->getMultigroup() == false){
			ng = 1;
		}

		int cw = mesh->getCellWidth();
		int ch = mesh->getCellHeight();
		size1 = cw*ch*ng;
		size2 = 4 + ng;
		cw = createAMPhi(size1, size2, ch*cw*ng);

	}
}

/**
 * cmfd Destructor clears all memory
 */
Cmfd::~Cmfd() {
}


/**
 * Create the loss matrix (A), fission matrix (A), and flux vector (phi)
 * @param number of columns needed in A matrix
 * @param number of rows needed in A matrix and M matrix
 * @param number of rows needed in phi vector
 */ 
int Cmfd::createAMPhi(PetscInt size1, PetscInt size2, int cells){

	int petsc_err = 0;

	petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &_A);
	size2 = size2 - 4;
	petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &_M);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, cells, &_phi_new);
	CHKERRQ(petsc_err);

	return petsc_err;
}

/** Computes the cross section for all MeshCells in the Mesh 
 * Create cross sections and fluxes for each cmfd cell by
 * energy condensing and volume averaging cross sections from
 * the MOC sweep.
 * @param pointer to an array of fsrs 
 *
 */
void Cmfd::computeXS(FlatSourceRegion* fsrs){

	_mesh->splitCorners();
	_flat_source_regions = fsrs;

	/* initialize variables */
	double volume, flux, abs, tot, nu_fis, chi;
	double* scat;
	double* mat_mult;
	double* mat_mult_a;
	double abs_tally_cell, nu_fis_tally_cell, dif_tally_cell, rxn_tally_cell;
	double vol_tally_cell, tot_tally_cell;
	double nu_fis_tally, dif_tally, rxn_tally, abs_tally, tot_tally;
	double scat_tally_cell[NUM_ENERGY_GROUPS];
	double scat_tot;

	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	Material* material;

	log_printf(DEBUG, "Enter Cmfd::computeXS(..)");

	/* For each mesh cell, we compute homogenized xs */
	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++)
	{
		meshCell = _mesh->getCells(i);

		/* Zeroes tallies for this energy group */
		abs_tally = 0.0;
		nu_fis_tally = 0.0;
		dif_tally = 0.0;
		tot_tally = 0.0;
		rxn_tally = 0.0;

		/* Computs flux weighted xs for each energy group */
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) 
		{
			abs_tally_cell = 0;
			nu_fis_tally_cell = 0;
			dif_tally_cell = 0;
			rxn_tally_cell = 0;
			vol_tally_cell = 0;
			tot_tally_cell = 0;
			scat_tot = 0;

			for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
			{
				scat_tally_cell[g] = 0;
			}

			/* loop over FSRs in mesh cell */
			std::vector<int>::iterator iter;
			for (iter = meshCell->getFSRs()->begin(); 
				 iter != meshCell->getFSRs()->end(); ++iter)
			{

				fsr = &_flat_source_regions[*iter];

				/* Gets FSR specific data. */
				material = fsr->getMaterial();
				chi = material->getChi()[e];
				volume = fsr->getVolume();
				//flux = fsr->getOldFlux()[e];
				flux = fsr->getFlux()[e];
				abs = material->getSigmaA()[e];
				tot = material->getSigmaT()[e];
				nu_fis = material->getNuSigmaF()[e];
				scat = material->getSigmaS();
				mat_mult = fsr->getMatMult();
				mat_mult_a = fsr->getMatMultA();

				abs_tally_cell += abs * flux * volume * mat_mult_a[e];
				tot_tally_cell += tot * flux * volume;
				//dif_tally_cell += flux  * volume / (3.0 * tot);
				nu_fis_tally_cell += nu_fis * flux * volume;
				rxn_tally_cell += flux * volume;
				vol_tally_cell += volume;


				for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
				{
					scat_tally_cell[g] += scat[g*NUM_ENERGY_GROUPS + e] 
						* flux * volume * mat_mult[g*NUM_ENERGY_GROUPS + e];
					scat_tot += scat[g * NUM_ENERGY_GROUPS + e];
					log_printf(DEBUG, "scattering from group %i to %i: %f", 
							   e, g, scat[g*NUM_ENERGY_GROUPS + e]);
				}

				if (chi >= meshCell->getChi()[e])
				{
					meshCell->setChi(chi,e);
				}

				if (nu_fis_tally_cell > 1e-10) /* Fission */
					dif_tally_cell += flux * volume / (3.0 * tot);
				else
					dif_tally_cell += flux * volume / (3.0 * tot - 2.0 * 
													   scat_tot);

			} /* end of looping through FSRs */

			/* For multi energy groups, we go ahead and set the xs for this 
			 * energy group */
			if (_mesh->getMultigroup() == true)
			{
				meshCell->setVolume(vol_tally_cell);
				meshCell->setSigmaA(abs_tally_cell / rxn_tally_cell, e);
				meshCell->setSigmaT(tot_tally_cell / rxn_tally_cell, e);
				meshCell->setNuSigmaF(nu_fis_tally_cell / rxn_tally_cell, e);
				meshCell->setDiffusivity(dif_tally_cell / rxn_tally_cell, e);
				meshCell->setOldFlux(rxn_tally_cell / vol_tally_cell, e);

				for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
				{
					meshCell->setSigmaS(scat_tally_cell[g] 
										/ rxn_tally_cell,e,g);
				}
			}
			/* For homogenized one energy group, we tally over all e's */
			else
			{
				abs_tally += abs_tally_cell;
				tot_tally += tot_tally_cell;
				nu_fis_tally += nu_fis_tally_cell;
				dif_tally += dif_tally_cell;
				rxn_tally += rxn_tally_cell;
			}
		}

		/* For homogenized one energy group, set xs after all e's are done */
		if (_mesh->getMultigroup() == false){
			meshCell->setVolume(vol_tally_cell);
			meshCell->setSigmaT(tot_tally / rxn_tally, 0);
			meshCell->setSigmaA(abs_tally / rxn_tally, 0);
			meshCell->setNuSigmaF(nu_fis_tally / rxn_tally, 0);
			meshCell->setDiffusivity(dif_tally / rxn_tally, 0);
			meshCell->setOldFlux(rxn_tally / vol_tally_cell, 0);
			meshCell->setChi(1,0);
		}
	}

	#if 0
	/* Computes keff = fission / abs based entirely on mat properties */
	double fis_tot = 0;
	double abs_tot = 0;
	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++){
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		{
			for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
			{
				/* FIXME: should there be a chi on the top of k? */
				fis_tot += meshCell->getChi()[e] * meshCell->getNuSigmaF()[g]
					* meshCell->getOldFlux()[g]*meshCell->getVolume();
			}
			abs_tot += meshCell->getSigmaA()[e]*meshCell->getOldFlux()[e]
				* meshCell->getVolume();
		}
	}

	/* print keff based on nu_fis / abs */
	log_printf(INFO, "fission rate / abs rate: %f", fis_tot / abs_tot);
	#endif
}

void Cmfd::computeDsxDirection(double x, double y, int e, MeshCell *meshCell, 
							   double d, double f, double flux, int cell_width)
{
	/* initialize variables */
	double d_next = 0, d_hat = 0, d_tilde = 0, 
		current = 0, flux_next = 0;
	MeshSurface *surf;
	MeshCell* meshCellNext;

	/* if cell on left side, set d_hat and d_tilde to 0 */
	if (x == 0)
	{
		if (_mesh->getBoundary(0) == REFLECTIVE)
		{
			meshCell->getMeshSurfaces(0)->setDDif(0.0, e);
		}
		else if (_mesh->getBoundary(0) == VACUUM)
		{
			surf = meshCell->getMeshSurfaces(0);
			current = - surf->getCurrent(e);

			/* set d_dif */
			surf->setDDif(2 * d / meshCell->getWidth() / 
						  (1 + 4 * d / meshCell->getWidth()), e);
			d_hat = 2 * d*f / meshCell->getWidth() 
				/ (1 + 4 * d*f / meshCell->getWidth());
			d_tilde = - (d_hat * flux + current 
						 / meshCell->getHeight()) / flux;
		}
	}
	/* if cell has a left neighbor, computes D's regularly */
	else
	{
		/* get mesh cell to left */
		meshCellNext = _mesh->getCells(y*cell_width + x - 1);
		d_next = meshCellNext->getDiffusivity()[e];
		flux_next = meshCellNext->getOldFlux()[e];

		/* set d_dif */
		MeshSurface *surf = meshCell->getMeshSurfaces(0);
		surf->setDDif(2.0 * d * d_next 
					  / (meshCell->getWidth() * d 
						 + meshCellNext->getWidth() * d_next), e);

		/* get diffusion correction term for meshCellNext */
		//f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

		/* compute d_hat */
		d_hat = 2.0 * d*f * d_next / 
			(meshCell->getWidth() * d * f 
			 + meshCellNext->getWidth() * d_next);

		/* Computes current: increment by outwards current on 
		 * next cell's RHS, decrement by outward current on LHS */
		current = 0.0;
		current += meshCellNext->getMeshSurfaces(2)->getCurrent(e);
		current -= meshCell->getMeshSurfaces(0)->getCurrent(e);

		/* compute d_tilde */
		d_tilde = -(d_hat * (flux - flux_next) + current  
					/ meshCell->getHeight()) / (flux_next + flux);

	}

	/* if abs(d_tilde) > abs(d_hat), make them equal in magnitude */
	if (fabs(d_tilde) > fabs(d_hat))
	{
		log_printf(DEBUG, "correcting Ds: LEFT group: %i, x: %i,"
				   " y: %i, dh: %f, dt: %f, c:%f", 
				   e, x, y, d_hat, d_tilde, current);

		/* d_tilde is positive */
		if (1 - fabs(d_tilde)/d_tilde < 1e-8)
		{
			d_hat   = - current/(2*flux*meshCell->getHeight());
			d_tilde = - current/(2*flux*meshCell->getHeight());
		}
		else
		{
			d_hat   = current/(2*flux_next*meshCell->getHeight());
			d_tilde = - current/(2*flux_next*meshCell->getHeight());
		}
	}


	log_printf(DEBUG, "cell: %i, group: %i, side: LEFT,"
			   " current: %f, dhat: %f, dtilde: %f", 
			   y*cell_width + x, e, current, d_hat, d_tilde);

	/* set d_hat and d_tilde */
	meshCell->getMeshSurfaces(0)->setDHat(d_hat, e);
	meshCell->getMeshSurfaces(0)->setDTilde(d_tilde, e);
}

/* compute the xs for all MeshCells in the Mesh */
void Cmfd::computeDs(){

	/* initialize variables */
	double d = 0, d_next = 0, d_hat = 0, d_tilde = 0, 
		current = 0, flux = 0, flux_next = 0, f = 1; //, f_next = 1;
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int ng = NUM_ENERGY_GROUPS;

	if (_mesh->getMultigroup() == false){
		ng = 1;
		_mesh->computeTotCurrents();
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();

	/* loop over all mesh cells */
	for (int y = 0; y < cell_height; y++)
	{
		for (int x = 0; x < cell_width; x++)
		{
			meshCell = _mesh->getCells(y*cell_width + x);

			for (int e = 0; e < ng; e++)
			{

				/* get diffusivity and flux for mesh cell */
				d = meshCell->getDiffusivity()[e];
				flux = meshCell->getOldFlux()[e];

				/* get diffusion correction term for meshCell */
				//f = computeDiffCorrect(d, meshCell->getWidth());
				f = 1.0;

				/* LEFT */
				computeDsxDirection(x, y, e, meshCell, d, f, flux, cell_width);

				/* RIGHT */
				/* if cell on right side, set d_hat and d_tilde to 0 */
				if (x == cell_width - 1){
					if (_mesh->getBoundary(2) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(2)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(2) == VACUUM){
						current = meshCell->getMeshSurfaces(2)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(2)->setDDif(2 * d / meshCell->getWidth() / (1 + 4 * d / meshCell->getWidth()), e);

						d_hat = 2 * d*f / meshCell->getWidth() / (1 + 4 * d*f / meshCell->getWidth());
						d_tilde = (d_hat * flux - current / meshCell->getHeight()) / flux;

					}
				}
				else{

					/* get mesh cell to the right */
					meshCellNext = _mesh->getCells(y*cell_width + x + 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(2)->setDDif(2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					//f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

					/* compute d_hat */
					d_hat = 2.0 * d * f * d_next / 
						(meshCell->getWidth() * d * f 
						 + meshCellNext->getWidth() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on right side */
					current += meshCell->getMeshSurfaces(2)->getCurrent(e);

					/* decrement current by outward current on next cell's left side */
					current -= meshCellNext->getMeshSurfaces(0)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current / meshCell->getHeight()) / (flux_next + flux);


				}

				log_printf(DEBUG, "cell: %i, group: %i, side:  RIGHT, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: RIGHT  group: %i, x: %i, y: %i, dh: %f, dt: %f, c: %f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux_next*meshCell->getHeight());
						d_tilde = - current/(2*flux_next*meshCell->getHeight());
					}
					else{
						d_hat   = current/(2*flux*meshCell->getHeight());
						d_tilde = - current/(2*flux*meshCell->getHeight());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(2)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(2)->setDTilde(d_tilde, e);

//////////////////////////////////////////////////////////////////////////////////////////////////



				/* BOTTOM */
				/* if cell on bottom side, set d_hat and d_tilde to 0 */

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getHeight());

				if (y == cell_height - 1){
					if (_mesh->getBoundary(1) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(1)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(1) == VACUUM){
						current = meshCell->getMeshSurfaces(1)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(1)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);

						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
						d_tilde = (d_hat * flux - current / meshCell->getWidth()) / flux;

					}
				}
				else{
					/* get mesh cell below */
					meshCellNext = _mesh->getCells((y+1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(1)->setDDif(2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					//f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next / 
						(meshCell->getHeight() * d * f 
						 + meshCellNext->getHeight() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on bottom side */
					current += meshCell->getMeshSurfaces(1)->getCurrent(e);

					/* decrement current by outward current on next cell's top side */
					current -= meshCellNext->getMeshSurfaces(3)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current / meshCell->getWidth()) / (flux_next + flux);

				}

				log_printf(DEBUG, "cell: %i, group: %i, side: BOTTOM, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: BOTTOM group: %i, x: %i, y: %i, dh: %f, dt: %f, c:%f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux_next*meshCell->getWidth());
						d_tilde = - current/(2*flux_next*meshCell->getWidth());
					}
					else{
						d_hat   = current/(2*flux*meshCell->getWidth());
						d_tilde = - current/(2*flux*meshCell->getWidth());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(1)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(1)->setDTilde(d_tilde, e);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////


				/* TOP */

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getHeight());

				/* if cell on top side, set d_hat and d_tilde to 0 */
				if (y == 0){
					if (_mesh->getBoundary(3) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(3)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(3) == VACUUM){
						current = - meshCell->getMeshSurfaces(3)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(3)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);

						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
						d_tilde = - (d_hat * flux + current / meshCell->getWidth()) / flux;
					}
				}
				else{
					/* get mesh cell above */
					meshCellNext = _mesh->getCells((y-1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(3)->setDDif(2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					//f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next / 
						(meshCell->getHeight() * d * f 
						 + meshCellNext->getHeight() * d_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on next cell's bottom side */
					current += meshCellNext->getMeshSurfaces(1)->getCurrent(e);

					/* decrement current by outward current on top side */
					current -= meshCell->getMeshSurfaces(3)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current / meshCell->getWidth()) / (flux_next + flux);
				}

				log_printf(DEBUG, "cell: %i, group: %i, side:    TOP, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: TOP    group: %i, x: %i, y: %i, dh: %f, dt: %f, c:%f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux*meshCell->getWidth());
						d_tilde = - current/(2*flux*meshCell->getWidth());
					}
					else{
						d_hat   = current/(2*flux_next*meshCell->getWidth());
						d_tilde = - current/(2*flux_next*meshCell->getWidth());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(3)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(3)->setDTilde(d_tilde, e);
			}
		}
	}
}

/* This is the original version of computes Ds */
void Cmfd::computeDsBackup(){

	/* initialize variables */
	double d = 0, d_next = 0, d_hat = 0, d_tilde = 0, 
		current = 0, flux = 0, flux_next = 0, f = 1, f_next = 1;
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int ng = NUM_ENERGY_GROUPS;

	if (_mesh->getMultigroup() == false){
		ng = 1;
		_mesh->computeTotCurrents();
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();

	/* loop over all mesh cells */
	for (int y = 0; y < cell_height; y++){
		for (int x = 0; x < cell_width; x++){
			meshCell = _mesh->getCells(y*cell_width + x);

			for (int e = 0; e < ng; e++){

				/* get diffusivity and flux for mesh cell */
				d = meshCell->getDiffusivity()[e];
				flux = meshCell->getOldFlux()[e];

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getWidth());

				/* LEFT */
				/* if cell on left side, set d_hat and d_tilde to 0 */
				if (x == 0)
				{
					if (_mesh->getBoundary(0) == REFLECTIVE)
					{
						meshCell->getMeshSurfaces(0)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(0) == VACUUM)
					{
						MeshSurface *surf = meshCell->getMeshSurfaces(0);
						current = - surf->getCurrent(e);

						/* set d_dif */
						surf->setDDif(2 * d / meshCell->getWidth() / 
									  (1 + 4 * d / meshCell->getWidth()), e);
						d_hat = 2 * d*f / meshCell->getWidth() 
							/ (1 + 4 * d*f / meshCell->getWidth());
						d_tilde = - (d_hat * flux + current 
									 / meshCell->getHeight()) / flux;
					}
				}
				/* if cell has a left neighbor, computes D's regularly */
				else
				{
					/* get mesh cell to left */
					meshCellNext = _mesh->getCells(y*cell_width + x - 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					MeshSurface *surf = meshCell->getMeshSurfaces(0);
					surf->setDDif(2.0 * d * d_next 
								  / (meshCell->getWidth() * d 
									 + meshCellNext->getWidth() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					f_next = computeDiffCorrect(d_next, 
												meshCellNext->getWidth());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next*f_next / 
						(meshCell->getWidth() * d * f 
						 + meshCellNext->getWidth() * d_next*f_next);

					/* Computes current: increment by outwards current on 
					 * next cell's RHS, decrement by outward current on LHS */
					current = 0.0;
					current += meshCellNext->getMeshSurfaces(2)->getCurrent(e);
					current -= meshCell->getMeshSurfaces(0)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current  
								/ meshCell->getHeight()) / (flux_next + flux);

				}

				/* if abs(d_tilde) > abs(d_hat), make them equal in magnitude */
				if (fabs(d_tilde) > fabs(d_hat))
				{
					log_printf(DEBUG, "correcting Ds: LEFT group: %i, x: %i,"
							   " y: %i, dh: %f, dt: %f, c:%f", 
							   e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8)
					{
						d_hat   = - current/(2*flux*meshCell->getHeight());
						d_tilde = - current/(2*flux*meshCell->getHeight());
					}
					else
					{
						d_hat   = current/(2*flux_next*meshCell->getHeight());
						d_tilde = - current/(2*flux_next*meshCell->getHeight());
					}
				}

				log_printf(DEBUG, "cell: %i, group: %i, side: LEFT,"
						   " current: %f, dhat: %f, dtilde: %f", 
						   y*cell_width + x, e, current, d_hat, d_tilde);

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(0)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(0)->setDTilde(d_tilde, e);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/* BOTTOM */
				/* if cell on bottom side, set d_hat and d_tilde to 0 */

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getHeight());

				if (y == cell_height - 1){
					if (_mesh->getBoundary(1) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(1)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(1) == VACUUM){
						current = meshCell->getMeshSurfaces(1)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(1)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);

						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
						d_tilde = (d_hat * flux - current / meshCell->getWidth()) / flux;

					}
				}
				else{
					/* get mesh cell below */
					meshCellNext = _mesh->getCells((y+1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(1)->setDDif(2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next*f_next / (meshCell->getHeight() * d*f + meshCellNext->getHeight() * d_next*f_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on bottom side */
					current += meshCell->getMeshSurfaces(1)->getCurrent(e);

					/* decrement current by outward current on next cell's top side */
					current -= meshCellNext->getMeshSurfaces(3)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current / meshCell->getWidth()) / (flux_next + flux);

				}

				log_printf(DEBUG, "cell: %i, group: %i, side: BOTTOM, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: BOTTOM group: %i, x: %i, y: %i, dh: %f, dt: %f, c:%f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux_next*meshCell->getWidth());
						d_tilde = - current/(2*flux_next*meshCell->getWidth());
					}
					else{
						d_hat   = current/(2*flux*meshCell->getWidth());
						d_tilde = - current/(2*flux*meshCell->getWidth());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(1)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(1)->setDTilde(d_tilde, e);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/* RIGHT */

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getWidth());

				/* if cell on right side, set d_hat and d_tilde to 0 */
				if (x == cell_width - 1){
					if (_mesh->getBoundary(2) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(2)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(2) == VACUUM){
						current = meshCell->getMeshSurfaces(2)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(2)->setDDif(2 * d / meshCell->getWidth() / (1 + 4 * d / meshCell->getWidth()), e);

						d_hat = 2 * d*f / meshCell->getWidth() / (1 + 4 * d*f / meshCell->getWidth());
						d_tilde = (d_hat * flux - current / meshCell->getHeight()) / flux;

					}
				}
				else{

					/* get mesh cell to the right */
					meshCellNext = _mesh->getCells(y*cell_width + x + 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(2)->setDDif(2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next*f_next / (meshCell->getWidth() * d*f + meshCellNext->getWidth() * d_next*f_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on right side */
					current += meshCell->getMeshSurfaces(2)->getCurrent(e);

					/* decrement current by outward current on next cell's left side */
					current -= meshCellNext->getMeshSurfaces(0)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux_next - flux) + current / meshCell->getHeight()) / (flux_next + flux);


				}

				log_printf(DEBUG, "cell: %i, group: %i, side:  RIGHT, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: RIGHT  group: %i, x: %i, y: %i, dh: %f, dt: %f, c: %f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux_next*meshCell->getHeight());
						d_tilde = - current/(2*flux_next*meshCell->getHeight());
					}
					else{
						d_hat   = current/(2*flux*meshCell->getHeight());
						d_tilde = - current/(2*flux*meshCell->getHeight());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(2)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(2)->setDTilde(d_tilde, e);

//////////////////////////////////////////////////////////////////////////////////////////////////

				/* TOP */

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getHeight());

				/* if cell on top side, set d_hat and d_tilde to 0 */
				if (y == 0){
					if (_mesh->getBoundary(3) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(3)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(3) == VACUUM){
						current = - meshCell->getMeshSurfaces(3)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(3)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);

						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
						d_tilde = - (d_hat * flux + current / meshCell->getWidth()) / flux;
					}
				}
				else{
					/* get mesh cell above */
					meshCellNext = _mesh->getCells((y-1)*cell_width + x);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(3)->setDDif(2.0 * d * d_next / (meshCell->getHeight() * d + meshCellNext->getHeight() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next*f_next / (meshCell->getHeight() * d*f + meshCellNext->getHeight() * d_next*f_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on next cell's bottom side */
					current += meshCellNext->getMeshSurfaces(1)->getCurrent(e);

					/* decrement current by outward current on top side */
					current -= meshCell->getMeshSurfaces(3)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current / meshCell->getWidth()) / (flux_next + flux);
				}

				log_printf(DEBUG, "cell: %i, group: %i, side:    TOP, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(INFO, "correcting Ds: TOP    group: %i, x: %i, y: %i, dh: %f, dt: %f, c:%f", e, x, y, d_hat, d_tilde, current);

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux*meshCell->getWidth());
						d_tilde = - current/(2*flux*meshCell->getWidth());
					}
					else{
						d_hat   = current/(2*flux_next*meshCell->getWidth());
						d_tilde = - current/(2*flux_next*meshCell->getWidth());
					}
				}

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(3)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(3)->setDTilde(d_tilde, e);
			}
		}
	}


	/* compute fis and abs rates based on tallied fluxes an XSs */
	double fis_tot = 0;
	double abs_tot = 0;
	double leak = 0;

	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++){
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < ng; e++){
			for (int g = 0; g < ng; g++){
				fis_tot += meshCell->getChi()[e]*meshCell->getNuSigmaF()[g]*meshCell->getOldFlux()[g]*meshCell->getVolume();
			}
			abs_tot += meshCell->getSigmaA()[e]*meshCell->getOldFlux()[e]*meshCell->getVolume();

			/* leakage */
			for (int s = 0; s < 4; s++){
				if (meshCell->getMeshSurfaces(s)->getBoundary() == VACUUM){
					if (s == 0){
						leak += meshCell->getMeshSurfaces(s)->getDHat()[e] * meshCell->getOldFlux()[e]*meshCell->getHeight();
						leak += meshCell->getMeshSurfaces(s)->getDTilde()[e] * meshCell->getOldFlux()[e]*meshCell->getHeight();
					}
					else if (s == 2){
						leak += meshCell->getMeshSurfaces(s)->getDHat()[e] * meshCell->getOldFlux()[e]*meshCell->getHeight();
						leak -= meshCell->getMeshSurfaces(s)->getDTilde()[e] * meshCell->getOldFlux()[e]*meshCell->getHeight();
					}
					else if (s == 1){
						leak += meshCell->getMeshSurfaces(s)->getDHat()[e] * meshCell->getOldFlux()[e]*meshCell->getWidth();
						leak -= meshCell->getMeshSurfaces(s)->getDTilde()[e] * meshCell->getOldFlux()[e]*meshCell->getWidth();
					}
					else if (s == 3){
						leak += meshCell->getMeshSurfaces(s)->getDHat()[e] * meshCell->getOldFlux()[e]*meshCell->getWidth();
						leak += meshCell->getMeshSurfaces(s)->getDTilde()[e] * meshCell->getOldFlux()[e]*meshCell->getWidth();
					}
				}
			}

		}
	}
}

/* Computes _quad_src based on (m+1/2) results */
void Cmfd::computeQuadSrc()
{
	/* Initializations */
	MeshSurface *s[4];
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int ng = NUM_ENERGY_GROUPS;
	double out[ng][8], in[ng][8];

	if (_mesh->getMultigroup() == false){
		ng = 1;
		//_mesh->computeTotCurrents();
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();

	/* loop over all mesh cells */
	for (int y = 0; y < cell_height; y++)
	{
		for (int x = 0; x < cell_width; x++)
		{
			meshCell = _mesh->getCells(y*cell_width + x);

			/* get four surfaces */
			for (int i = 0; i < 4; i++) 
				s[i] = meshCell->getMeshSurfaces(i);

			for (int e = 0; e < ng; e++)
			{
				out[e][0] = s[2]->getFlux(e,0);
				out[e][1] = s[1]->getFlux(e,0);
				out[e][2] = s[3]->getFlux(e,1);
				out[e][3] = s[2]->getFlux(e,1);
				out[e][4] = s[0]->getFlux(e,0);
				out[e][5] = s[3]->getFlux(e,0);
				out[e][6] = s[1]->getFlux(e,1);
				out[e][7] = s[0]->getFlux(e,1);
			}


			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[0].flux[0] = %f",
					   x, y, s[0]->getFlux(0, 0));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[0].flux[1] = %f",
					   x,y, s[0]->getFlux(0, 1));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[2].flux[0] = %f",
					   x,y, s[2]->getFlux(0, 0));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[2].flux[1] = %f",
					   x,y, s[2]->getFlux(0, 1));

			if (x == 0)
			{
				if (_mesh->getBoundary(0) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][5] = s[0]->getFlux(e,1);
						in[e][6] = s[0]->getFlux(e,0);
					}
				}
				else if (_mesh->getBoundary(0) == VACUUM)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][5] = 0;
						in[e][6] = 0;
					}
				}
			}
			else
			{
				meshCellNext = _mesh->getCells(y*cell_width + x - 1);
				for (int e = 0; e < ng; e++)
				{
					in[e][5] = meshCellNext->getMeshSurfaces(2)->getFlux(e,0);
					in[e][6] = meshCellNext->getMeshSurfaces(2)->getFlux(e,1);
				}			
			}


			if (x == cell_width - 1)
			{
				if (_mesh->getBoundary(2) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][1] = s[2]->getFlux(e,1);
						in[e][2] = s[2]->getFlux(e,0);
					}
				}
				else if (_mesh->getBoundary(2) == VACUUM)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][1] = 0;
						in[e][2] = 0;
					}
				}
			}
			else
			{
				meshCellNext = _mesh->getCells(y*cell_width + x + 1);
				for (int e = 0; e < ng; e++)
				{
					in[e][1] = meshCellNext->getMeshSurfaces(0)->getFlux(e,0);
					in[e][2] = meshCellNext->getMeshSurfaces(0)->getFlux(e,1);
				}			
			}			
			
			if (y == 0)
			{
				if (_mesh->getBoundary(3) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][3] = s[3]->getFlux(e,0);
						in[e][4] = s[3]->getFlux(e,1);
					}
				}
				else if (_mesh->getBoundary(3) == VACUUM)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][3] = 0;
						in[e][4] = 0;
					}
				}
			}
			else
			{
				meshCellNext = _mesh->getCells( (y - 1) * cell_width + x);
				for (int e = 0; e < ng; e++)
				{
					in[e][3] = meshCellNext->getMeshSurfaces(1)->getFlux(e,1);
					in[e][4] = meshCellNext->getMeshSurfaces(1)->getFlux(e,0);
				}			
			}

			if (y == cell_height - 1)
			{
				if (_mesh->getBoundary(1) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][7] = s[1]->getFlux(e,0);
						in[e][0] = s[1]->getFlux(e,1);
					}
				}
				else if (_mesh->getBoundary(1) == VACUUM)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][7] = 0;
						in[e][0] = 0;
					}
				}
			}
			else
			{
				meshCellNext = _mesh->getCells( (y + 1) * cell_width + x);
				for (int e = 0; e < ng; e++)
				{
					in[e][7] = meshCellNext->getMeshSurfaces(3)->getFlux(e,1);
					in[e][0] = meshCellNext->getMeshSurfaces(0)->getFlux(e,0);
				}			
			}

			/* Now that we have all the in's and out's, computes src */
			double l = meshCell->getL();
			for (int e = 0; e < ng; e++)
			{
				double xs = meshCell->getSigmaT()[e];
				double ex = exp(-xs * l);
				double sum_quad_flux = 0;

				for (int i = 0; i < 8; i++)
				{
					double src = xs * (out[e][i] - ex * in[e][i]) / (1.0 - ex);
					meshCell->setQuadSrc(src, e, i);

					sum_quad_flux += src/xs - (out[e][i] - in[e][i])/(xs * l);
				}
				meshCell->setSumQuadFlux(sum_quad_flux, e);
			}
		}
	}
}	 

/*
* CMFD solver that solves the diffusion problem
* @param solve methed - either diffusion or cmfd (acceleration)
* @param iteration number of in MOC solver - used for plotting
* @return k-effective
*/
double Cmfd::computeCMFDFluxPower(solveType solveMethod, int moc_iter){

	log_printf(INFO, "Running diffusion solver...");

	/* Initializes variables */
	MeshCell* meshCell;
	int ng = NUM_ENERGY_GROUPS;
	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int petsc_err;
	Vec phi_old;
	PetscInt size1;
	PetscScalar sumold, sumnew, scale_val, eps;
	PetscReal rtol = 1e-10;
	PetscReal atol = 1e-10;
	double criteria = FLUX_CONVERGENCE_THRESH;
	Vec sold, snew, res;
	int max_outer = 100;


	/* If running single group acceleration, set number of group to 1 */
	if (_mesh->getMultigroup() == false){
		ng = 1;
	}

	/* Constructs array to store old flux */
	size1 = cw*ch*ng;
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, size1, &phi_old);
	CHKERRQ(petsc_err);

	/* construct A matrix, M matrix, and flux vector */
	petsc_err = MatZeroEntries(_A);
	petsc_err = MatZeroEntries(_M);
	petsc_err = constructAMPhi(_A, _M, phi_old, solveMethod);
	CHKERRQ(petsc_err);

	if (solveMethod == DIFFUSION){
		max_outer = 1000;
	}


	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &sold);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &snew);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &res);
	CHKERRQ(petsc_err);

	/* Assemblies vectors and matrices */
	petsc_err = VecAssemblyBegin(phi_old);
	petsc_err = VecAssemblyEnd(phi_old);
	petsc_err = VecAssemblyBegin(_phi_new);
	petsc_err = VecAssemblyEnd(_phi_new);
	petsc_err = VecAssemblyBegin(sold);
	petsc_err = VecAssemblyEnd(sold);
	petsc_err = VecAssemblyBegin(snew);
	petsc_err = VecAssemblyEnd(snew);
	petsc_err = VecAssemblyBegin(res);
	petsc_err = VecAssemblyEnd(res);
	CHKERRQ(petsc_err);
	petsc_err = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyBegin(_M, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(_M, MAT_FINAL_ASSEMBLY);
	CHKERRQ(petsc_err);

	/* Sets _phi_new to phi_old */
	petsc_err = VecCopy(phi_old, _phi_new);
	CHKERRQ(petsc_err);

	/* Initializes KSP solver */
	KSP ksp;
	petsc_err = KSPCreate(PETSC_COMM_WORLD, &ksp);
	petsc_err = KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, PETSC_DEFAULT);
	petsc_err = KSPSetType(ksp, KSPGMRES);
	petsc_err = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	petsc_err = KSPSetOperators(ksp, _A, _A, SAME_NONZERO_PATTERN);
	petsc_err = KSPSetUp(ksp);
	petsc_err = KSPSetFromOptions(ksp);
	CHKERRQ(petsc_err);

	/* Obtains initial source and finds intitial k_eff */
	petsc_err = MatMult(_M, _phi_new, snew);
	petsc_err = VecSum(snew, &sumnew);
	petsc_err = MatMult(_A, _phi_new, sold);
	petsc_err = VecSum(sold, &sumold);
	_keff = float(sumnew)/float(sumold);
	log_printf(INFO, "CMFD iter: %i, keff: %f", 0, _keff);
	CHKERRQ(petsc_err);

	/* Normalizes initial source */
	petsc_err = MatMult(_M, _phi_new, sold);
	petsc_err = VecSum(sold, &sumold);
	scale_val = (cw * ch * ng) / sumold;
	petsc_err = VecScale(sold, scale_val);
	CHKERRQ(petsc_err);

	sumold = cw * ch * ng;
	int iter = 0;

	for (iter = 0; iter < max_outer; iter++){

		petsc_err = KSPSolve(ksp, sold, _phi_new);
		petsc_err = MatMult(_M, _phi_new, snew);
		petsc_err = VecSum(snew, &sumnew);
		CHKERRQ(petsc_err);

		/* Computes the new keff */
		_keff = sumnew / sumold;
		log_printf(INFO, "CMFD iter: %i, keff: %f", iter + 1, _keff);
		petsc_err = VecScale(sold, _keff);

		/* Computes the L2 norm of source error: eps */
		scale_val = 1e-15;
		petsc_err = VecShift(snew, scale_val);
		petsc_err = VecPointwiseDivide(res, sold, snew);
		scale_val = -1;
		petsc_err = VecShift(res, scale_val);
		CHKERRQ(petsc_err);
		petsc_err = VecNorm(res, NORM_2, &eps);

		/* Scales the fissino source such that vector sum is \# cells x ng */
		eps = eps / (cw * ch * ng);
		scale_val = (cw * ch * ng) / sumnew;
		petsc_err = VecScale(snew, scale_val);
		CHKERRQ(petsc_err);

		petsc_err = VecCopy(snew, sold);
		CHKERRQ(petsc_err);

		if (iter > 5 && eps < criteria){
			break;
		}
	}

	/* Rescales the new and old flux */
	double factor = 1;
	petsc_err = MatMult(_M, phi_old, sold);
	petsc_err = VecSum(sold, &sumold);
	scale_val = factor * (cw * ch * ng) / sumold;
	petsc_err = VecScale(phi_old, scale_val);
	CHKERRQ(petsc_err);

	petsc_err = MatMult(_M, _phi_new, snew);
	petsc_err = VecSum(snew, &sumnew);
	scale_val = factor * (cw * ch * ng) / sumnew;
	petsc_err = VecScale(_phi_new, scale_val);
	CHKERRQ(petsc_err);


	/* Passes new and old flux into Mesh Cell to set flux */
	PetscScalar *old_phi;
	PetscScalar *new_phi;
	petsc_err = VecGetArray(phi_old, &old_phi);
	petsc_err = VecGetArray(_phi_new, &new_phi);
	CHKERRQ(petsc_err);
	for (int i = 0; i < cw*ch; i++)
	{
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < ng; e++)
		{
			meshCell->setOldFlux(double(old_phi[i*ng + e]), e);
			meshCell->setNewFlux(double(new_phi[i*ng + e]), e);
		}
		
		log_printf(INFO, "When CMFD converges, old phi = %f, new phi = %f", 
				   double(old_phi[i*ng+ng-1]), double(new_phi[i*ng+ng-1]));
	}

	petsc_err = VecRestoreArray(phi_old, &old_phi);
	petsc_err = VecRestoreArray(_phi_new, &new_phi);
	CHKERRQ(petsc_err);

	/*
	MatMult(_A,_phi_new, sold);
	sumold = 1/_keff;
	MatScale(_M,sumold);
	MatMult(_M,_phi_new,snew);
	sumold = -1;
	VecWAXPY(res, sumold, snew, sold);
	VecSum(res,&sumold);

	log_printf(DEBUG, "CMFD/DIFFUSION residual: %f", double(sumold));
	*/

	/* destroy matrices and vectors */
	petsc_err = VecDestroy(&phi_old);
	petsc_err = VecDestroy(&snew);
	petsc_err = VecDestroy(&sold);
	petsc_err = VecDestroy(&res);
	petsc_err = KSPDestroy(&ksp);
	CHKERRQ(petsc_err);

	std::string string;

	if (solveMethod == DIFFUSION){
		if (_plotter->plotDiffusion() == true){
			string = "diff";
			_plotter->plotCMFDflux(_mesh, string, moc_iter);
		}
	}

	if (_plotter->plotKeff()){
		_plotter->plotCMFDKeff(_mesh, moc_iter);
	}

	if (_plotter->plotCurrent()){
		string = "cmfd";
		_plotter->plotCMFDflux(_mesh, string, moc_iter);
	}

	if (solveMethod == CMFD){
		_mesh->setKeffCMFD(_keff, moc_iter);
	}
	return _keff;
}

/* Computes the flux in each mesh cell using LOO */
double Cmfd::computeLooFluxPower(solveType solveMethod, int moc_iter){
	log_printf(INFO, "Running low order MOC solver...");
	int iter, max_outer = 1000; 
	if (solveMethod == DIFFUSION){
		max_outer = 1000;
	}
	double old_keff = 1.0;

	/* Obtains info about the meshes */
	MeshCell* meshCell;
	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int ng = NUM_ENERGY_GROUPS;
	if (_mesh->getMultigroup() == false){
		ng = 1;
	}

	/* Copies old fluxes into new fluxes */
	for (int i = 0; i < cw * ch; i++)
	{
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < ng; e++) 
			meshCell->setNewFlux(meshCell->getOldFlux()[e], e);
	}

	/* Pre-compute factors that do not change between iterations */
	double **sum_quad_flux, **quad_xs, **ratio, **expo, **tau;
	sum_quad_flux = new double*[cw*ch];
	quad_xs = new double*[cw*ch];
	ratio = new double*[cw*ch];
	expo = new double*[cw*ch];
	tau = new double*[cw*ch];

	double l = _mesh->getCells(0)->getL();

	for (int i = 0; i < cw * ch; i++)
	{
		sum_quad_flux[i] = new double[NUM_ENERGY_GROUPS];
		quad_xs[i] = new double[NUM_ENERGY_GROUPS];
		ratio[i] = new double[NUM_ENERGY_GROUPS];
		expo[i] = new double[NUM_ENERGY_GROUPS];
		tau[i] = new double[NUM_ENERGY_GROUPS];
		
		for (int e = 0; e < ng; e++)
		{
			sum_quad_flux[i][e] = 0;
			quad_xs[i][e] = _mesh->getCells(i)->getSigmaT()[e];
			tau[i][e] = quad_xs[i][e] * l;
			expo[i][e] = exp(-tau[i][e]);
			ratio[i][e] = (1 - expo[i][e]) / tau[i][e];
		}
	}

	_keff = 1.0;

	/* Starts LOO acceleration iteration, we do not update src, quad_src, 
	 * quad_flux, old_flux, as they are computed from the MOC step (i.e., 
	 * order m+1/2) and should not be updated during acceleration step. */
	for (iter = 0; iter < max_outer; iter++)
	{
		old_keff = _keff;

		/* Computes new cell averaged source, looping over energy groups */
		double **new_src;
		new_src = new double*[cw*ch];
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			new_src[i] = new double[NUM_ENERGY_GROUPS];
			for (int e = 0; e < ng; e++)
			{
				new_src[i][e] = 0.0;
				for (int g = 0; g < ng; g++)
				{
					new_src[i][e] += meshCell->getSigmaS()[e*ng+g];
					new_src[i][e] += meshCell->getChi()[e] 
						* meshCell->getNuSigmaF()[g] / _keff;
				}
				new_src[i][e] *= meshCell->getNewFlux()[e];
				log_printf(DEBUG, "Cell averaged source for cell %d, energy %d"
						   " is %e", i, e, new_src[i][e]);
			}
		}

		/* Updates 8 quadrature sources based on form function */
		double **new_quad_src;
		new_quad_src = new double*[cw*ch];
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			new_quad_src[i] = new double[8 * NUM_ENERGY_GROUPS];
			for (int e = 0; e < ng; e++)
			{
				for (int g = 0; g < 8; g++)
				{
					int d = e * ng + g;
					/* Notice getSrc()[e] returns the $\bar{Q}_g^{(m)}$ */
					new_quad_src[i][d] = meshCell->getQuadSrc()[d] 
						* new_src[i][e] / meshCell->getSrc()[e];	
					log_printf(DEBUG, "Old Mesh Averaged Source = %e", 
							   meshCell->getSrc()[e]);
					log_printf(DEBUG, "Updated quadrature source for cell %d" 
							   " energy %d, track %d is %f", i, e, g, 
							   new_quad_src[i][d]);
				}
			}
		}

		/* Sweeps over geometry, solve LOO MOC */
		/* FIXME: we do not really accumulate scalar flux as we do in the outter
		 * MOC sweeps right? We only deal with angular flux in this case? */
		for (int e = 0; e < ng; e++)
		{
			double flux;
			int i, g, d;
			int i_array[]  = {2,3,1,1,1,0,2,2};
			int g_array[]  = {0,5,0,2,4,1,4,6};
			int i_array2[] = {3,3,1,0,0,0,2,3};
			int g_array2[] = {0,2,7,2,4,6,3,6};

			/* 1st loop */
			flux = _mesh->getCells(2)->getQuadFlux()[e*ng + 0];		
			for (int x = 0; x < 8; x++)
			{
				for (int y = 0; y < 8; y++)
				{
					i = i_array[x];
					g = g_array[y];
					d = e * ng + g;
					sum_quad_flux[i][e] += flux * ratio[i][e] + 
						new_quad_src[i][d] * l * ratio[i][e];
					flux = expo[i][e] * flux + new_quad_src[i][d] * l 
						* ratio[i][e];
				}
			}
			_mesh->getCells(2)->setQuadFlux(flux, e, 0);
			flux = _mesh->getCells(2)->getQuadFlux()[e*ng + 7];		
			for (int x = 7; x >=0 ; x--)
			{
				for (int y = 7; y >= 0; y--)
				{
					i = i_array[x];
					g = g_array[y];
					d = e * ng + g;
					sum_quad_flux[i][e] += flux * ratio[i][e] + 
						new_quad_src[i][d] * l * ratio[i][e];
					flux = expo[i][e] * flux + new_quad_src[i][d] * l 
						* ratio[i][e];
				}
			}
			_mesh->getCells(2)->setQuadFlux(flux, e, 7);


			/* 2nd loop */
			flux = _mesh->getCells(3)->getQuadFlux()[e*ng + 0];		
			for (int x = 0; x < 8; x++)
			{
				for (int y = 0; y < 8; y++)
				{
					i = i_array2[x];
					g = g_array2[y];
					d = e * ng + g;
					sum_quad_flux[i][e] += flux * ratio[i][e] + 
						new_quad_src[i][d] * l * ratio[i][e];
					flux = expo[i][e] * flux + new_quad_src[i][d] * l 
						* ratio[i][e];
				}
			}
			_mesh->getCells(3)->setQuadFlux(flux, e, 0);
			flux = _mesh->getCells(3)->getQuadFlux()[e*ng + 7];		
			for (int x = 7; x >=0 ; x--)
			{
				for (int y = 7; y >= 0; y--)
				{
					i = i_array2[x];
					g = g_array2[y];
					d = e * ng + g;
					sum_quad_flux[i][e] += flux * ratio[i][e] + 
						new_quad_src[i][d] * l * ratio[i][e];
					flux = expo[i][e] * flux + new_quad_src[i][d] * l 
						* ratio[i][e];
				}
			}
			_mesh->getCells(3)->setQuadFlux(flux, e, 7);
		}		

		/* Computes new cell-averaged scalar flux based on new_sum_quad_flux */
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			for (int e = 0; e < ng; e++) 
			{
				meshCell->setNewFlux(meshCell->getOldFlux()[e] 
									 * sum_quad_flux[i][e] 
									 / meshCell->getSumQuadFlux()[e], e);
			}
		}

		/* Checks for convergence */
		/* Computes keff assuming zero leakage */
		double fis_tot = 0, abs_tot = 0;
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			for (int e = 0; e < ng; e++)
			{
				for (int g = 0; g < ng; g++)
				{
					/* FIXME: should there be a chi on the top of k? */
					fis_tot += meshCell->getChi()[e]
						* meshCell->getNuSigmaF()[g]
						* meshCell->getNewFlux()[g] * meshCell->getVolume();
				}
				abs_tot += meshCell->getSigmaA()[e] * meshCell->getNewFlux()[e]
					* meshCell->getVolume();
			}
		}
		_keff = fis_tot / abs_tot; 
		log_printf(NORMAL, "%d-th LOO iteration k = %f", iter, _keff);

		if ((iter > 5) && (fabs(_keff - old_keff) / _keff < 1e-6))
		{
			/* new flux is already in place */
			/* just need to print stuff */
			std::string string;
			if (solveMethod == DIFFUSION)
			{
				if (_plotter->plotDiffusion() == true){
					string = "diff";
					_plotter->plotCMFDflux(_mesh, string, moc_iter);
				}
			}

			if (_plotter->plotKeff())
			{
				_plotter->plotCMFDKeff(_mesh, moc_iter);
			}

			if (_plotter->plotCurrent())
			{
				string = "loo";
				_plotter->plotCMFDflux(_mesh, string, moc_iter);
			}

			if (solveMethod == CMFD)
			{
				_mesh->setKeffCMFD(_keff, moc_iter);
			}

			return _keff;
		}
	}
	
	log_printf(WARNING, "Keff not converging after %d iterations", max_outer);
	return 1.0;
}


int Cmfd::constructAMPhi(Mat A, Mat M, Vec phi_old, solveType solveMethod){

	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int ng;
	MeshCell* meshCell;
	int petsc_err = 0;

	PetscInt indice1, indice2;
	PetscScalar value;

	ng = NUM_ENERGY_GROUPS;
	if (_mesh->getMultigroup() == false){
		ng = 1;
	}

	/* loop over mesh cells in y direction */
	for (int y = 0; y < ch; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cw; x++){

			/* loop over energy groups */
			for (int e = 0; e < ng; e++){

				/* get mesh cell */
				meshCell = _mesh->getCells(y*cw + x);

				/* get old flux */
				indice1 = int((y*cw + x)*ng+e);
				petsc_err = VecSetValues(phi_old, 1, &indice1, &meshCell->getOldFlux()[e], INSERT_VALUES);

				CHKERRQ(petsc_err);

				/* diagonal - A: add outscattering term to diagonal */
				value = meshCell->getSigmaA()[e] * meshCell->getVolume();
				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* add out-scattering term to diagonal */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = meshCell->getSigmaS()[e*ng + g] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng+e;
						petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
						CHKERRQ(petsc_err);
					}
				}


				/* diagonal - M: add fission terms to diagonal in M */
				for (int g = 0; g < ng; g++){
					value = meshCell->getChi()[e] * meshCell->getNuSigmaF()[g] * meshCell->getVolume();
					indice1 = (y*cw + x)*ng+e;
					indice2 = (y*cw + x)*ng + g;
					petsc_err = MatSetValues(M, 1, &indice1, 1, &indice2, &value, INSERT_VALUES);
				CHKERRQ(petsc_err);
				}


				/* scattering in */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = - meshCell->getSigmaS()[g*ng + e] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng + g;
						petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
						CHKERRQ(petsc_err);
					}
				}

				/* RIGHT */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(2)->getDHat()[e]      - meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(2)->getDDif()[e] * meshCell->getHeight();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1,1 , &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* transport in */
				if (x != cw - 1){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(2)->getDHat()[e] + meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(2)->getDDif()[e] * meshCell->getHeight();

					indice1 = (y*cw + x)*ng + e;
					indice2 = (y*cw + x + 1)*ng + e;
					petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
					CHKERRQ(petsc_err);
				}

				/* LEFT */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(0)->getDHat()[e]      + meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(0)->getDDif()[e] * meshCell->getHeight();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* transport in */
				if (x != 0){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(0)->getDHat()[e] - meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(0)->getDDif()[e] * meshCell->getHeight();

					indice1 = (y*cw + x)*ng + e;
					indice2 = (y*cw + x - 1)*ng + e;
					petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
					CHKERRQ(petsc_err);
				}

				/* BELOW */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(1)->getDHat()[e]      - meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(1)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* transport in */
				if (y != ch - 1){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(1)->getDHat()[e] + meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(1)->getDDif()[e] * meshCell->getWidth();

					indice1 = (y*cw + x)*ng + e;
					indice2 = ((y+1)*cw + x)*ng + e;
					petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
					CHKERRQ(petsc_err);
				}

				/* ABOVE */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(3)->getDHat()[e]      + meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 =  (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* transport in */
				if (y != 0){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(3)->getDHat()[e] - meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

					indice1 = (y*cw + x)*ng + e;
					indice2 = ((y-1)*cw + x)*ng + e;
					petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, 
											 &value, ADD_VALUES);
					CHKERRQ(petsc_err);
				}
			}
		}
	}

	return petsc_err;
}

/* Updates the MOC flux in each FSR based on results from the acceleration
 * steps; we use the ratio of the cell averaged flux to update the FSR flux
 * @param MOC iteration number
 */
void Cmfd::updateMOCFlux(int iteration){

	log_printf(INFO, "Updating MOC flux using form func from acceleration...");

	/* initialize variables */
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	double old_flux, new_flux;//, fsr_new_flux;
	double* flux;

	int cw = _mesh->getCellWidth();
	int ch = _mesh->getCellHeight();
	int ng = NUM_ENERGY_GROUPS;

	/* loop over mesh cells */
	for (int i = 0; i < cw * ch; i++){

		/* get pointer to current mesh cell */
		meshCell = _mesh->getCells(i);

		for (int e = 0; e < ng; e++){

			if (_mesh->getMultigroup() == true){
				old_flux = meshCell->getOldFlux()[e];
				new_flux = meshCell->getNewFlux()[e];
			}
			else {
				old_flux = meshCell->getOldFlux()[0];
				new_flux = meshCell->getNewFlux()[0];
			}

			/* Report new flux and old flux   */
			if (e == NUM_ENERGY_GROUPS - 1)
			{
				log_printf(NORMAL, "Update flux in Cell: %i,"
						   " old = %f, new = %f, new/old = %f", 
						   i, old_flux, new_flux, new_flux / old_flux);
			}

			/* loop over FRSs in mesh cell */
			std::vector<int>::iterator iter;
			for (iter = meshCell->getFSRs()->begin(); 
				 iter != meshCell->getFSRs()->end(); ++iter) 
			{
				fsr = &_flat_source_regions[*iter];

				/* get fsr flux */
				flux = fsr->getFlux();

				flux[e] = new_flux / old_flux * flux[e];
			}
		}
	}
}



double Cmfd::computeDiffCorrect(double d, double h){

	/* compute correction - F */
	/* double alpha, mu, expon;
	double rho, F;
	rho = 0.0;
	for (int p = 0; p < NUM_POLAR_ANGLES; p++){
		mu = std::cos(std::asin(_quad->getSinTheta(p)));
		expon = exp(- h / (3 * d * mu));
		alpha = (1 + expon) / (1 - expon) - 2 * mu / h;
		rho += mu * _quad->getWeight(p) * alpha;
	}

    F = 1 + h * rho / (2 * d);
	*/

	return 1.0;

}

Mat Cmfd::getA(){
	return _A;
}

Mat Cmfd::getM(){
	return _M;
}

Vec Cmfd::getPhiNew(){
	return _phi_new;

}

double Cmfd::getKeff(){
	return _keff;
}


/* Store the mesh cell averaged source before a MOC sweep */
void Cmfd::storePreMOCMeshSource(FlatSourceRegion* fsrs)
{
	_flat_source_regions = fsrs;

	/* initialize variables */
	double volume, source;
	double vol_tally_cell;
	double source_tally_cell, source_tally;

	MeshCell* meshCell;
	FlatSourceRegion* fsr;

	log_printf(DEBUG, "Enter Cmfd::storePreMOCMeshSource(..)");

	/* For each mesh cell, we compute homogenized xs */
	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++)
	{
		meshCell = _mesh->getCells(i);

		/* Zeroes tallies for this energy group */
		source_tally = 0;

		/* Computes flux weighted xs for each energy group */
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) 
		{
			vol_tally_cell = 0;
			source_tally_cell = 0;

			/* loop over FSRs in mesh cell, accumulate cell tally */
			std::vector<int>::iterator iter;
			for (iter = meshCell->getFSRs()->begin(); 
				 iter != meshCell->getFSRs()->end(); ++iter)
			{
				fsr = &_flat_source_regions[*iter];
				volume = fsr->getVolume();
				source = fsr->getSource()[e];

				source_tally_cell += source * volume;
				vol_tally_cell += volume;
			} 

			/* For multi energy groups, we go ahead and set the xs for this 
			 * energy group */
			if (_mesh->getMultigroup() == true)
				meshCell->setSrc(source_tally_cell / vol_tally_cell, e);
			else /* For homogenized one energy group, we tally over all e's */
				source_tally += source_tally_cell;
		}

		/* For homogenized one energy group, set xs after all e's are done */
		if (_mesh->getMultigroup() == false)
			meshCell->setSrc(source_tally / vol_tally_cell, 0);
	}
}

/* Set the old flux for each FSR equal to FSR flux */
void Cmfd::setOldFSRFlux()
{
	/* initialize variables */
	FlatSourceRegion* fsr;
	int num_fsrs = _geom->getNumFSRs();

	/* Compute total fission source for this region */
	for (int r = 0; r < num_fsrs; r++) 
	{
		/* Get pointers to important data structures */
		fsr = &_flat_source_regions[r];
	  
		/* loop over energy groups */
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		{
			fsr->setOldFlux(e, fsr->getFlux()[e]);
		  
		}
	}
}


/* set pointer to array of fsrs
 * @param pointer to arrary of fsrs
 */
void Cmfd::setFSRs(FlatSourceRegion* fsrs)
{
  _flat_source_regions = fsrs;
} 

