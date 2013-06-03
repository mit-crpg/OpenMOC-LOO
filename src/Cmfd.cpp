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
Cmfd::Cmfd(Geometry* geom, Plotter* plotter, Mesh* mesh, bool runCmfd, 
		   bool useDiffusionCorrection, double l2_norm_conv_thresh,
		   TrackGenerator *track_generator) {

	_geom = geom;
	_plotter = plotter;
	_mesh = mesh;
	_quad = new Quadrature(TABUCHI);

	_num_azim = track_generator->getNumAzim();
	_spacing = track_generator->getSpacing();
	_l2_norm = 1.0;
	_l2_norm_conv_thresh = l2_norm_conv_thresh;
	_use_diffusion_correction = useDiffusionCorrection;

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
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, cells, &_source_old);
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
void Cmfd::computeXS(){

	/* split corner currents to side surfaces */
	/* FIXME: implement splitcorners for quad currents for LOO */
	_mesh->splitCorners();

	/* initialize variables */
	double volume, flux, abs, tot, nu_fis, chi;
	double* scat;
	double abs_tally_group, nu_fis_tally_group, dif_tally_group, 
		rxn_tally_group, vol_tally_group, tot_tally_group;
	double nu_fis_tally = 0, dif_tally = 0, rxn_tally = 0, abs_tally = 0, 
		tot_tally = 0, scat_tally = 0;
	double scat_tally_group[NUM_ENERGY_GROUPS];

	std::vector<int>::iterator iter;
	int i,e,g;

	/* create pointers to objects */
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	Material* material;

	/* loop over mesh cells */
#if USE_OPENMP
#pragma omp parallel for private(i, e, g, volume, flux, abs, tot, nu_fis, \
    chi, scat, abs_tally_group, nu_fis_tally_group, \
    dif_tally_group, rxn_tally_group, vol_tally_group, tot_tally_group, \
    nu_fis_tally, dif_tally, rxn_tally, abs_tally, tot_tally, \
    scat_tally_group, iter, meshCell, fsr, material)
#endif 

	for (i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++){
		meshCell = _mesh->getCells(i);

		/* if single group, zero tallies */
		if (_mesh->getMultigroup() == false)
		{
			abs_tally = 0.0;
			nu_fis_tally = 0.0;
			dif_tally = 0.0;
			tot_tally = 0.0;
			rxn_tally = 0.0;
			scat_tally = 0.0;
		}

		/* loop over energy groups */
		for (e = 0; e < NUM_ENERGY_GROUPS; e++) 
		{

			/* zero tallies for this group */
			abs_tally_group = 0;
			nu_fis_tally_group = 0;
			dif_tally_group = 0;
			rxn_tally_group = 0;
			vol_tally_group = 0;
			tot_tally_group = 0;

			/* zero each group to group scattering tally */
			for (g = 0; g < NUM_ENERGY_GROUPS; g++)
			{
				scat_tally_group[g] = 0;
			}

			/* loop over FSRs in mesh cell */
			for (iter = meshCell->getFSRs()->begin(); 
				 iter != meshCell->getFSRs()->end(); ++iter)
			{
				fsr = &_flat_source_regions[*iter];

				/* Gets FSR specific data. */
				volume = fsr->getVolume();
				flux = fsr->getFlux()[e];
				material = fsr->getMaterial();
				chi = material->getChi()[e];
				abs = material->getSigmaA()[e];
				tot = material->getSigmaT()[e];
				nu_fis = material->getNuSigmaF()[e];
				scat = material->getSigmaS();


				/* increment tallies for this group */
				abs_tally_group += abs * flux * volume;
				tot_tally_group += tot * flux * volume;
				nu_fis_tally_group += nu_fis * flux * volume;
				rxn_tally_group += flux * volume;
				vol_tally_group += volume;
				dif_tally_group += flux * volume / (3.0 * tot);

				/* increment group to group scattering tallies */
				for (g = 0; g < NUM_ENERGY_GROUPS; g++){
					/* scattering from group e into g */
					scat_tally_group[g] += scat[ g * NUM_ENERGY_GROUPS + e] 
						* flux * volume;
				}

				/* choose a chi for this group */
				if (chi >= meshCell->getChi()[e]){
					meshCell->setChi(chi,e);
				}
			}

			/* if multigroup, set the multigroup parameters */
			if (_mesh->getMultigroup()){
				meshCell->setVolume(vol_tally_group);
				meshCell->setSigmaA(abs_tally_group / rxn_tally_group, e);
				meshCell->setSigmaT(tot_tally_group / rxn_tally_group, e);
				meshCell->setNuSigmaF(nu_fis_tally_group / rxn_tally_group, e);
				meshCell->setDiffusivity(dif_tally_group / rxn_tally_group, e);
				meshCell->setOldFlux(rxn_tally_group / vol_tally_group, e);

				for (g = 0; g < NUM_ENERGY_GROUPS; g++)
				{
					/* This means SigmaS[e * ng + g] = $\Sigma_{s,e\to g}$. */
					meshCell->setSigmaS(scat_tally_group[g] / rxn_tally_group,
										e, g);
				}
			}
			/* if single group, add group-wise tallies up */
			else{
				abs_tally += abs_tally_group;
				tot_tally += tot_tally_group;
				nu_fis_tally += nu_fis_tally_group;
				dif_tally += dif_tally_group;
				rxn_tally += rxn_tally_group;
				for (g = 0; g < NUM_ENERGY_GROUPS; g++)
					scat_tally += scat_tally_group[g];
			}
		}

		/* if single group, set single group parameters */
		if (_mesh->getMultigroup() == false){
			meshCell->setVolume(vol_tally_group);
			meshCell->setSigmaT(tot_tally / rxn_tally, 0);
			meshCell->setSigmaA(abs_tally / rxn_tally, 0);
			meshCell->setNuSigmaF(nu_fis_tally / rxn_tally, 0);
			meshCell->setDiffusivity(dif_tally / rxn_tally, 0);
			meshCell->setOldFlux(rxn_tally / vol_tally_group, 0);
			log_printf(INFO, "mesh = %d, rxn tally = %.10f, vol = %.10f,"
					   " mesh's old flux = %.10f", 
					   i, rxn_tally, vol_tally_group, 
					   meshCell->getOldFlux()[0]);
			meshCell->setChi(1, 0);
			/* FIXME: SG needs to add up all scattering */
			meshCell->setSigmaS(scat_tally / rxn_tally, 0, 0);
		}
	}
}

void Cmfd::computeDsxDirection(double x, double y, int e, MeshCell *meshCell, 
							   double d, double f, double flux, int cell_width, 
							   double dt_weight)
{
	/* initialize variables */
	double d_next = 0, d_hat = 0, d_tilde = 0, 
		current = 0, flux_next = 0, f_next = 0;
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
		f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

		/* compute d_hat */
		d_hat = 2.0 * d*f * d_next * f_next / 
			(meshCell->getWidth() * d * f 
			 + meshCellNext->getWidth() * d_next * f_next);

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
		log_printf(DEBUG, "correcting Ds: LEFT group: %i, x: %f,"
				   " y: %f, dh: %f, dt: %f, c:%f", 
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


	log_printf(DEBUG, "cell: %f, group: %i, side: LEFT,"
			   " current: %.10f, dhat: %f, dtilde: %f", 
			   y*cell_width + x, e, current, d_hat, d_tilde);

	/* set d_hat and d_tilde */
	d_tilde = meshCell->getMeshSurfaces(0)->getDTilde()[e] 
		* (1 - dt_weight) + dt_weight * d_tilde;
	meshCell->getMeshSurfaces(0)->setDHat(d_hat, e);
	meshCell->getMeshSurfaces(0)->setDTilde(d_tilde, e);
}

/* compute the xs for all MeshCells in the Mesh */
void Cmfd::computeDs(){

	/* initialize variables */
	double d = 0, d_next = 0, d_hat = 0, d_tilde = 0, 
		current = 0, flux = 0, flux_next = 0, f = 1, f_next = 1;
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int ng = NUM_ENERGY_GROUPS;
	int x, y, e;

	if (_mesh->getMultigroup() == false){
		ng = 1;
		_mesh->computeTotCurrents();
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();
	double dt_weight = 0.66;

	/* loop over all mesh cells */
#if USE_OPENMP
#pragma omp parallel for private(meshCell, x, y, e, \
    d, flux, f, d_hat, d_tilde, current, d_next, \
    flux_next, meshCellNext)
#endif 
	for (y = 0; y < cell_height; y++)
	{
		for (x = 0; x < cell_width; x++)
		{
			meshCell = _mesh->getCells(y*cell_width + x);

			for (e = 0; e < ng; e++)
			{

				/* get diffusivity and flux for mesh cell */
				d = meshCell->getDiffusivity()[e];
				flux = meshCell->getOldFlux()[e];

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getWidth());

				/* LEFT */
				computeDsxDirection(x, y, e, meshCell, d, f, flux, cell_width,
									dt_weight);

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
					f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

					/* compute d_hat */
					d_hat = 2.0 * d * f * d_next * f_next / 
						(meshCell->getWidth() * d * f 
						 + meshCellNext->getWidth() * d_next * f_next);

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
				d_tilde = meshCell->getMeshSurfaces(2)->getDTilde()[e] 
					* (1 - dt_weight) + dt_weight * d_tilde;
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
					f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d * f * d_next * f_next / 
						(meshCell->getHeight() * d * f 
						 + meshCellNext->getHeight() * d_next * f_next);

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
				d_tilde = meshCell->getMeshSurfaces(1)->getDTilde()[e] 
					* (1 - dt_weight) + dt_weight * d_tilde;
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
					f_next = computeDiffCorrect(d_next, meshCellNext->getHeight());

					/* compute d_hat */
					d_hat = 2.0 * d * f * d_next * f_next / 
						(meshCell->getHeight() * d * f 
						 + meshCellNext->getHeight() * d_next * f_next);

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
				d_tilde = meshCell->getMeshSurfaces(3)->getDTilde()[e] 
					* (1 - dt_weight) + dt_weight * d_tilde;
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
					d_hat = 2.0 * d * f * d_next * f_next / 
						(meshCell->getWidth() * d * f 
						 + meshCellNext->getWidth() * d_next * f_next);

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
					d_hat = 2.0 * d*f * d_next* f_next 
						/ (meshCell->getHeight() * d * f + 
						   meshCellNext->getHeight() * d_next * f_next);

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
					d_hat = 2.0 * d*f * d_next * f_next 
						/ (meshCell->getWidth() * d * f 
						   + meshCellNext->getWidth() * d_next * f_next);

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
					d_hat = 2.0 * d*f * d_next * f_next 
						/ (meshCell->getHeight() * d * f 
						   + meshCellNext->getHeight() * d_next * f_next);

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

/* Computes _quad_flux based on _quad_current */
void Cmfd::computeQuadFlux()
{
	/* Initializations */
	MeshSurface *s[4];
	MeshCell* meshCell;
	int ng = NUM_ENERGY_GROUPS;

	if (_mesh->getMultigroup() == false){
		ng = 1;
		_mesh->computeTotQuadCurrents();
		// FIXME: debug
		_mesh->computeTotCurrents();
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();

	/* the factor that we devide everyone by is cos(45degree) * surface len */
	/* May need to fixme */
	double scale = _mesh->getCells(0)->getWidth() * SIN_THETA_45;

	/* loop over all mesh cells */
	for (int y = 0; y < cell_height; y++)
	{
		for (int x = 0; x < cell_width; x++)
		{
			meshCell = _mesh->getCells(y*cell_width + x);

			/* get four surfaces */
			for (int i = 0; i < 4; i++) 
			{
				s[i] = meshCell->getMeshSurfaces(i);
				for (int e = 0; e < ng; e++)
				{
					for (int ind = 0; ind < 2; ind ++)
					{
						s[i]->setQuadFlux(s[i]->getQuadCurrent(e, ind) / scale, 
										  e, ind);
					}
					/* For debugging purpose, compute partial currents */
					double current = s[i]->getQuadCurrent(e, 0) 
						+ s[i]->getQuadCurrent(e, 1);
					log_printf(ACTIVE, "cmfd generates %.10f,"
							   " LOO generates %.10f",
							   s[i]->getCurrent(e), current);


					s[i]->setCurrent(current, e);

					/* Prints to screen quad current and quad flux */
					log_printf(ACTIVE, "cell %d surface %d energy %d's current"
							   " is %.10f, quad fluxes = %.10f %.10f", 
							   y*cell_width+x, i, e, current,
							   s[i]->getQuadFlux(e, 0), s[i]->getQuadFlux(e,1));
				}
			}
		}
	}

	/* Debugging */
	log_printf(DEBUG, "set cell 2 side 1 track 0: %f", 
			   _mesh->getCells(2)->getMeshSurfaces(1)->getQuadFlux(0, 0));

	return;
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

	if (_mesh->getMultigroup() == false)
		ng = 1;

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
				out[e][0] = s[2]->getQuadFlux(e,0);
				out[e][1] = s[1]->getQuadFlux(e,0);
				out[e][2] = s[3]->getQuadFlux(e,1);
				out[e][3] = s[2]->getQuadFlux(e,1);
				out[e][4] = s[0]->getQuadFlux(e,0);
				out[e][5] = s[3]->getQuadFlux(e,0);
				out[e][6] = s[1]->getQuadFlux(e,1);
				out[e][7] = s[0]->getQuadFlux(e,1);
			}


			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[0].flux[0] = %f",
					   x, y, s[0]->getQuadFlux(0, 0));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[0].flux[1] = %f",
					   x,y, s[0]->getQuadFlux(0, 1));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[2].flux[0] = %f",
					   x,y, s[2]->getQuadFlux(0, 0));
			log_printf(DEBUG, "Cell (x,y) = (%d, %d), surface[2].flux[1] = %f",
					   x,y, s[2]->getQuadFlux(0, 1));

			if (x == 0)
			{
				if (_mesh->getBoundary(0) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][5] = s[0]->getQuadFlux(e,1);
						in[e][6] = s[0]->getQuadFlux(e,0);
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
					in[e][5] = meshCellNext->getMeshSurfaces(2)->getQuadFlux(e,0);
					in[e][6] = meshCellNext->getMeshSurfaces(2)->getQuadFlux(e,1);
				}			
			}


			if (x == cell_width - 1)
			{
				if (_mesh->getBoundary(2) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][1] = s[2]->getQuadFlux(e,1);
						in[e][2] = s[2]->getQuadFlux(e,0);
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
					in[e][1] = meshCellNext->getMeshSurfaces(0)->getQuadFlux(e,0);
					in[e][2] = meshCellNext->getMeshSurfaces(0)->getQuadFlux(e,1);
				}			
			}			
			
			if (y == 0)
			{
				if (_mesh->getBoundary(3) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][3] = s[3]->getQuadFlux(e,0);
						in[e][4] = s[3]->getQuadFlux(e,1);
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
					in[e][3] = meshCellNext->getMeshSurfaces(1)->getQuadFlux(e,1);
					in[e][4] = meshCellNext->getMeshSurfaces(1)->getQuadFlux(e,0);
				}			
			}

			if (y == cell_height - 1)
			{
				if (_mesh->getBoundary(1) == REFLECTIVE)
				{
					for (int e = 0; e < ng; e++)
					{
						in[e][7] = s[1]->getQuadFlux(e,0);
						in[e][0] = s[1]->getQuadFlux(e,1);
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
					in[e][7] = meshCellNext->getMeshSurfaces(3)->getQuadFlux(e,1);
					in[e][0] = meshCellNext->getMeshSurfaces(0)->getQuadFlux(e,0);
				}			
			}

			/* Now that we have all the in's and out's, computes src */
			double l = meshCell->getL();
			for (int e = 0; e < ng; e++)
			{
				double xs = meshCell->getSigmaT()[e];
				double ex = exp(- xs * l);
				double sum_quad_flux = 0;

				for (int t = 0; t < 8; t++)
				{
					double src = xs * (out[e][t] - ex * in[e][t]) / (1.0 - ex);
					meshCell->setQuadSrc(src, e, t);

					sum_quad_flux += src/xs + (in[e][t] - out[e][t])/(xs * l);
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
double Cmfd::computeCMFDFluxPower(solveType solveMethod, int moc_iter)
{

	log_printf(INFO, "Running diffusion solver...");

	/* initialize variables */
	MeshCell* meshCell;
	int ng = NUM_ENERGY_GROUPS;
	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int petsc_err;
	Vec phi_old, sold, snew, res;
	PetscInt size;
	int max_outer, iter = 0;
	PetscScalar sumold, sumnew, scale_val, eps;
	PetscReal rtol = 1e-10;
	PetscReal atol = 1e-10;
	std::string string;

	/* if single group, set ng (number of groups) to 1 */
	if (_mesh->getMultigroup() == false){
		ng = 1;
	}

	/* create petsc array to store old flux */
	petsc_err = VecCreate(PETSC_COMM_WORLD, &phi_old);
	size = ch * cw * ng;
	petsc_err = VecSetSizes(phi_old, PETSC_DECIDE, size);
	petsc_err = VecSetFromOptions(phi_old);
	CHKERRQ(petsc_err);

	/* zero out and construct memory efficient versions of
	 * A matrix, M matrix, and flux vector */
	petsc_err = MatZeroEntries(_A);
	petsc_err = MatZeroEntries(_M);
	petsc_err = constructAMPhi(_A, _M, phi_old, solveMethod);
	CHKERRQ(petsc_err);

	/* if solve method is DIFFUSION, set max_outer to large number
	 * to allow solver to fully converge */
	if (solveMethod == DIFFUSION){
		max_outer = 1000;
	}
	/* if solve method is CMFD, set max_outer such that flux
	 * partially converges */
	else{
		max_outer = 100;
	}

	/* create old source and residual vectors */
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &sold);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &snew);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &res);
	CHKERRQ(petsc_err);

	/* assembly vectors and matrices*/
	petsc_err = VecAssemblyBegin(phi_old);
	petsc_err = VecAssemblyEnd(phi_old);
	petsc_err = VecAssemblyBegin(_phi_new);
	petsc_err = VecAssemblyEnd(_phi_new);
	petsc_err = VecAssemblyBegin(_source_old);
	petsc_err = VecAssemblyEnd(_source_old);
	petsc_err = VecAssemblyBegin(sold);
	petsc_err = VecAssemblyEnd(sold);
	petsc_err = VecAssemblyBegin(snew);
	petsc_err = VecAssemblyEnd(snew);
	petsc_err = VecAssemblyBegin(res);
	petsc_err = VecAssemblyEnd(res);
	petsc_err = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyBegin(_M, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(_M, MAT_FINAL_ASSEMBLY);
	CHKERRQ(petsc_err);

	/* set _phi_new to phi_old */
	petsc_err = VecCopy(phi_old, _phi_new);
	CHKERRQ(petsc_err);

	/* initialize KSP solver */
	KSP ksp;
	petsc_err = KSPCreate(PETSC_COMM_WORLD, &ksp);
	petsc_err = KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, PETSC_DEFAULT);
	petsc_err = KSPSetType(ksp, KSPGMRES);
	petsc_err = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	petsc_err = KSPSetOperators(ksp, _A, _A, SAME_NONZERO_PATTERN);
	petsc_err = KSPSetUp(ksp);
	petsc_err = KSPSetFromOptions(ksp);
	CHKERRQ(petsc_err);

	/* get initial source and find initial k_eff */
	petsc_err = MatMult(_M, _phi_new, snew);
	petsc_err = VecSum(snew, &sumnew);
	petsc_err = MatMult(_A, _phi_new, sold);
	petsc_err = VecSum(sold, &sumold);
	_keff = float(sumnew)/float(sumold);
	log_printf(INFO, "CMFD iter: %i, keff: %f", iter, _keff);
	CHKERRQ(petsc_err);

	/* recompute and normalize initial source */
	petsc_err = MatMult(_M, _phi_new, sold);
	petsc_err = VecSum(sold, &sumold);
	scale_val = (cw * ch * ng) / sumold;
	petsc_err = VecScale(sold, scale_val);
	sumold = cw * ch * ng;
	CHKERRQ(petsc_err);

	/* diffusion solver */
	for (iter = 0; iter < max_outer; iter++){

		/* Solve x = A_inverse * b problem and compute new source */
		petsc_err = KSPSolve(ksp, sold, _phi_new);
		petsc_err = MatMult(_M, _phi_new, snew);
		petsc_err = VecSum(snew, &sumnew);
		CHKERRQ(petsc_err);

		/* compute and set keff */
		_keff = sumnew / sumold;
		log_printf(INFO, "CMFD iter: %i, keff: %f", iter + 1, _keff);
		petsc_err = VecScale(sold, _keff);

		/* compute the L2 norm of source error */
		scale_val = 1e-15;
		petsc_err = VecShift(snew, scale_val);
		petsc_err = VecPointwiseDivide(res, sold, snew);
		scale_val = -1;
		petsc_err = VecShift(res, scale_val);
		CHKERRQ(petsc_err);
		petsc_err = VecNorm(res, NORM_2, &eps);

		/* compute error */
		eps = eps / (cw * ch * ng);
		scale_val = (cw * ch * ng) / sumnew;
		petsc_err = VecScale(snew, scale_val);
		CHKERRQ(petsc_err);

		/* set old source to new source */
		petsc_err = VecCopy(snew, sold);
		CHKERRQ(petsc_err);

		/* check for convergence for the CMFD iterative solver using 
		 * _l2_norm_conv_thresh as criteria */
		if (iter > 5 && eps < _l2_norm_conv_thresh)
		{
			log_printf(DEBUG, " CMFD converges in %d iterations", iter);
			break;
		}
	}

	/* rescale the new and old flux */
	petsc_err = MatMult(_M, _phi_new, snew);
	petsc_err = VecSum(snew, &sumnew);
	scale_val = (cw * ch * ng) / sumnew;
	petsc_err = VecScale(_phi_new, scale_val);
	petsc_err = MatMult(_M, phi_old, sold);
	petsc_err = VecSum(sold, &sumold);
	scale_val = (cw * ch * ng) / sumold;
	petsc_err = VecScale(phi_old, scale_val);
	CHKERRQ(petsc_err);

	log_printf(INFO, " CMFD total fission source = %f", sumnew);


	PetscScalar *old_phi;
	PetscScalar *new_phi;
	petsc_err = VecGetArray(phi_old, &old_phi);
	petsc_err = VecGetArray(_phi_new, &new_phi);
	CHKERRQ(petsc_err);

	for (int i = 0; i < cw*ch; i++){
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < ng; e++){
			meshCell->setOldFlux(double(old_phi[i*ng + e]), e);
			meshCell->setNewFlux(double(new_phi[i*ng + e]), e);
		}
	}

	petsc_err = VecRestoreArray(phi_old, &old_phi);
	petsc_err = VecRestoreArray(_phi_new, &new_phi);
	CHKERRQ(petsc_err);

	/* Computes L2 norm between the source that enters the CMFD acceleration 
	 * step and the one coming out of converged CMFD step to decided whether 
	 * the outter MOC iteration / source iteration should quit */
	if (moc_iter > 0)
	{
		petsc_err = fisSourceNorm(snew, moc_iter);
	}

    /* Copies source new to source old */
	petsc_err = VecCopy(snew, _source_old);
	CHKERRQ(petsc_err);

	/* compute the residual */
	MatMult(_A,_phi_new, sold);
	sumold = 1/_keff;
	MatScale(_M,sumold);
	MatMult(_M,_phi_new,snew);
	sumold = -1;
	VecWAXPY(res, sumold, snew, sold);
	VecSum(res,&sumold);
	log_printf(DEBUG, "CMFD/DIFFUSION residual: %f", double(sumold));

	/* destroy matrices and vectors */
	petsc_err = VecDestroy(&phi_old);
	petsc_err = VecDestroy(&snew);
	petsc_err = VecDestroy(&sold);
	petsc_err = VecDestroy(&res);
	petsc_err = KSPDestroy(&ksp);
	CHKERRQ(petsc_err);


	/* plot flux, current, and k_eff */
	if (solveMethod == DIFFUSION){
		if (_plotter->plotDiffusion() == true){
			string = "diff";
			_plotter->plotCMFDflux(_mesh, string, moc_iter);
		}
	}

	if (_plotter->plotKeff()){
		_plotter->plotCMFDKeff(_mesh, moc_iter);
	}

	/* comments out plotting current for each iteration */
	/*
	if (_plotter->plotCurrent()){
		string = "cmfd";
		_plotter->plotCMFDflux(_mesh, string, moc_iter);
	}
	*/

	if (solveMethod == CMFD){
		_mesh->setKeffCMFD(_keff, moc_iter);
	}

	return _keff;
}

/* Computes the flux in each mesh cell using LOO */
double Cmfd::computeLooFluxPower(solveType solveMethod, int moc_iter, 
								 double k_MOC)
{
	log_printf(INFO, "Running low order MOC solver...");

	int iter, max_outer = 10; 

	if (solveMethod == DIFFUSION){
		max_outer = 1000;
	}

	_keff = k_MOC;

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


	double l = _mesh->getCells(0)->getL();

	/* Allocate memories for terms internal to the LOO iterative solver */
	double **sum_quad_flux, **quad_xs, **ratio, **expo, **tau, **new_src;
	double **new_quad_src;
	try
	{
		sum_quad_flux = new double*[cw*ch];
		quad_xs = new double*[cw*ch];
		ratio = new double*[cw*ch];
		expo = new double*[cw*ch];
		tau = new double*[cw*ch];
		new_src = new double*[cw*ch];
		new_quad_src = new double*[cw*ch];

		for (int i = 0; i < cw * ch; i++)
		{
			sum_quad_flux[i] = new double[NUM_ENERGY_GROUPS];
			quad_xs[i] = new double[NUM_ENERGY_GROUPS];
			ratio[i] = new double[NUM_ENERGY_GROUPS];
			expo[i] = new double[NUM_ENERGY_GROUPS];
			tau[i] = new double[NUM_ENERGY_GROUPS];
			new_src[i] = new double[NUM_ENERGY_GROUPS];
			new_quad_src[i] = new double[8 * NUM_ENERGY_GROUPS];
		}
	}
	catch (std::exception &e)
	{
		log_error("Unable to allocate memory for variables inside loo. "
				   "Backtrace: \n%s", e.what());
	}

	for (int i = 0; i < cw * ch; i++)
	{
		for (int e = 0; e < ng; e++)
		{
			new_src[i][e] = 0.0;
			sum_quad_flux[i][e] = 0;
			quad_xs[i][e] = _mesh->getCells(i)->getSigmaT()[e];
			/* Pre-computes cross-section related terms because they do 
			   not change between iterations */
			tau[i][e] = quad_xs[i][e] * l;
			expo[i][e] = exp(-tau[i][e]);
			ratio[i][e] = (1 - expo[i][e]) / tau[i][e];
		}
		for (int t = 0; t < 8 * ng; t++)
			new_quad_src[i][t] = 0.0;
	}

	double eps = 0.0;
	double old_flux[cw*ch], new_flux[cw*ch], xs;
	for (int i = 0; i < cw * ch; i++)
	{
		meshCell = _mesh->getCells(i);
		old_flux[i] = 0;
		/* integrates over energy */
		for (int e = 0; e < ng; e++)
		{
			xs = meshCell->getNuSigmaF()[e];
			old_flux[i] += xs * meshCell->getOldFlux()[e];
		} 
	}



	/* Starts LOO acceleration iteration, we do not update src, quad_src, 
	 * quad_flux, old_flux, as they are computed from the MOC step (i.e., 
	 * order m+1/2) and should not be updated during acceleration step. */
	for (iter = 0; iter < max_outer; iter++)
	{
		//old_keff = _keff;
		
		/* Resets terms to zeros for each LOO iteration */
		for (int i = 0; i < cw * ch; i++)
		{
			for (int e = 0; e < ng; e++)
			{
				sum_quad_flux[i][e] = 0.0;
				new_src[i][e] = 0.0;
				for (int t = 0; t < 8; t++)
				{
					int d = e * 8 + t;
					new_quad_src[i][d] = 0.0;
				}
			}
		}
		
		/* Computes new cell averaged source, looping over energy groups */
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);

			for (int e = 0; e < ng; e++)
			{
				for (int g = 0; g < ng; g++)
				{
					new_src[i][e] += meshCell->getSigmaS()[g*ng+e] 
						* meshCell->getNewFlux()[g] * ONE_OVER_FOUR_PI ;
					new_src[i][e] += meshCell->getChi()[e] *
						meshCell->getNuSigmaF()[g] / _keff 
						* meshCell->getNewFlux()[g] * ONE_OVER_FOUR_PI ;
				}
			}


			/*
			double fission_source = 0, scatter_source = 0;
			double *scalar_flux = meshCell->getNewFlux();
			double *nu_sigma_f = meshCell->getNuSigmaF();
			double *sigma_s = meshCell->getSigmaS();
			double *chi = meshCell->getChi();

			for (int e = 0; e < ng; e++)
				fission_source += nu_sigma_f[e] * scalar_flux[e];

			for (int e = 0; e < ng; e++)
			{
				scatter_source = 0;

				if (_mesh->getMultigroup())
				{
					for (int g = 0; g < ng; g++)
						scatter_source += sigma_s[g*ng+e] * scalar_flux[g];		
				}
	   
				new_src[i][e] = (fission_source * chi[e] / _keff 
								 + scatter_source) * ONE_OVER_FOUR_PI;
 
				if (e == ng - 1)
					log_printf(INFO, "Cell averaged source for cell %d,"
							   " energy %d is %e", i, e, new_src[i][e]);
			}
			*/


		} /* finishing looping over i; exit to iter level */

		double src_ratio;
		/* Updates 8 quadrature sources based on form function */
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			for (int e = 0; e < ng; e++)
			{
				/* getSrc()[e] returns the $\bar{Q}_g^{(m)}$ */
				src_ratio = new_src[i][e] / meshCell->getSrc()[e];
				log_printf(DEBUG, "Old Mesh Averaged Source = %e", 
						   meshCell->getSrc()[e]);
				/* FIXME */
				//src_ratio = 1.0;
				for (int t = 0; t < 8; t++)
				{
					int d = e * 8 + t;
					new_quad_src[i][d] = meshCell->getQuadSrc()[d] * src_ratio;
				}
				log_printf(ACTIVE, "Average source ratio for cell %d" 
						   " energy %d, by %.10f", i, e, src_ratio);
			}
		} /* finish iterating over i; exit to iter level */

		/* Sweeps over geometry, solve LOO MOC */
		for (int e = 0; e < ng; e++)
		{
			double flux;
			int i, t, d;
			/* Forward direction, i_array[] contains cell #, g_array[] contains
			 * track # */
			int i_array[]  = {2,3,1,1,1,0,2,2};
			//int i_array[]  = {3,2,0,0,0,1,3,3};
			int t_array[]  = {0,5,0,2,4,1,4,6};
		
			int i_array2[] = {3,3,1,0,0,0,2,3};
			//int i_array2[] = {2,2,0,1,1,1,3,2};
			int t_array2[] = {0,2,7,2,4,6,3,6};

			/* 1st loop */
			/* Get the initial angular flux */
			flux = _mesh->getCells(i_array[0])->getMeshSurfaces(1)->getQuadFlux(e,0);	
			log_printf(ACTIVE, " Forward, begin with %f", flux);
			for (int x = 0; x < 8; x++)
			{
				i = i_array[x];
				t = t_array[x];
				d = e * 8 + t;
				/* Accumulate angular flux to $\bar{\psi}_g^{8,(n+1)}$ */
				sum_quad_flux[i][e] += flux * ratio[i][e] + 
					new_quad_src[i][d] * (1- ratio[i][e]) / quad_xs[i][e];
				/* Update angular flux: $\psi_{out} = \psi_{in} *
				 * e^{-\Sigma L} + Q/\Sigma (1 - e^{-\Sigma L}) */
				flux -= (flux * tau[i][e] - new_quad_src[i][d] * l) 
					* ratio[i][e];
			}
			_mesh->getCells(i_array[0])->getMeshSurfaces(1)->setQuadFlux(flux, e, 0);
			log_printf(ACTIVE, "  end with %f", flux);

			flux = _mesh->getCells(i_array[0])->getMeshSurfaces(1)->getQuadFlux(e,1); 
			log_printf(ACTIVE, " Backward, begin with %f", flux);
			for (int x = 7; x >=0 ; x--)
			{
				i = i_array[x];
				t = t_array[x];
				d = e * 8 + t;
				/* Accumulate angular flux to $\bar{\psi}_g^{8,(n+1)}$ */
				sum_quad_flux[i][e] += flux * ratio[i][e] + 
					new_quad_src[i][d] * (1- ratio[i][e]) / quad_xs[i][e];
				/* Update angular flux: $\psi_{out} = \psi_{in} *
				 * e^{-\Sigma L} + Q/\Sigma (1 - e^{-\Sigma L}) */
				flux -= (flux * tau[i][e] - new_quad_src[i][d] * l) 
					* ratio[i][e];
			}
			_mesh->getCells(i_array[0])->getMeshSurfaces(1)->setQuadFlux(flux, e, 1);
			log_printf(ACTIVE, "  end with %f", flux);



			/* 2nd loop */
			flux = _mesh->getCells(i_array2[0])->getMeshSurfaces(1)->getQuadFlux(e, 0);
			log_printf(ACTIVE, " Forward, loop 2, begin with %f", flux);

			for (int x = 0; x < 8; x++)
			{
				i = i_array2[x];
				t = t_array2[x];
				d = e * 8 + t;
				/* Accumulate angular flux to $\bar{\psi}_g^{8,(n+1)}$ */
				sum_quad_flux[i][e] += flux * ratio[i][e] + 
					new_quad_src[i][d] * (1- ratio[i][e]) / quad_xs[i][e];
				/* Update angular flux: $\psi_{out} = \psi_{in} *
				 * e^{-\Sigma L} + Q/\Sigma (1 - e^{-\Sigma L}) */
				flux -= (flux * tau[i][e] - new_quad_src[i][d] * l) 
					* ratio[i][e];
			}
			_mesh->getCells(i_array2[0])->getMeshSurfaces(1)->setQuadFlux(flux, e, 0);
			log_printf(ACTIVE, "  end with %f", flux);

			flux = _mesh->getCells(i_array2[0])->getMeshSurfaces(1)->getQuadFlux(e, 1);
			log_printf(ACTIVE, " Backward, loop 2, begin with %f", flux);
			for (int x = 7; x >=0 ; x--)
			{
				i = i_array2[x];
				t = t_array2[x];
				d = e * 8 + t;
				/* Accumulate angular flux to $\bar{\psi}_g^{8,(n+1)}$ */
				sum_quad_flux[i][e] += flux * ratio[i][e] + 
					new_quad_src[i][d] * (1 - ratio[i][e]) / quad_xs[i][e];
				/* Update angular flux: $\psi_{out} = \psi_{in} *
				 * e^{-\Sigma L} + Q/\Sigma (1 - e^{-\Sigma L}) */
				flux -= (flux * tau[i][e] - new_quad_src[i][d] * l) 
					* ratio[i][e];
			}
			_mesh->getCells(i_array2[0])->getMeshSurfaces(1)->setQuadFlux(flux, e, 1);
			log_printf(ACTIVE, "  end with %f", flux);
		} /* finish looping over energy; exit to iter level */				

		double phi_ratio;
		/* Computes new cell-averaged scalar flux based on new_sum_quad_flux */
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			for (int e = 0; e < ng; e++) 
			{
				/* FIXME: 1/8 is the weight on angular quadrature */
				phi_ratio =  sum_quad_flux[i][e] 
					/ meshCell->getSumQuadFlux()[e]; //  / 8;

				/* FIXME */
				//phi_ratio = 1.0;

				meshCell->setNewFlux(meshCell->getOldFlux()[e] * phi_ratio, e);

				log_printf(ACTIVE, "Update the cell-averaged scalar flux by "
						   "factor of %.10f, old phi = %.10f", 
						   phi_ratio, meshCell->getOldFlux()[e]);
			}
		}

		/* Computes keff assuming zero leakage */
		double fis_tot = 0, abs_tot = 0, vol_tot = 0;
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			vol_tot += meshCell->getVolume();
			for (int e = 0; e < ng; e++)
			{
				for (int g = 0; g < ng; g++)
				{
					fis_tot += meshCell->getChi()[e]
						* meshCell->getNuSigmaF()[g]
						* meshCell->getNewFlux()[g] ;//* meshCell->getVolume();
				}
		    
				abs_tot += meshCell->getSigmaA()[e] * meshCell->getNewFlux()[e]
					;//* meshCell->getVolume();
			}
		}
		_keff = fis_tot / abs_tot; 

		double normalize_factor = 1.0 / fis_tot * ch * cw;
		/* Normalizes flux based on fission source FIXME: update? */
		for (int i = 0; i < cw * ch; i++)
		{
			meshCell = _mesh->getCells(i);
			for (int e = 0; e < ng; e++)
			{
				meshCell->setNewFlux(meshCell->getNewFlux()[e] 
									 * normalize_factor, e);
			}
		}
		fis_tot *= normalize_factor; 
		abs_tot *= normalize_factor;

		log_printf(ACTIVE, "Normalize all scalar flux by %f", normalize_factor);

		/* Computes the L2 norm of point-wise-division of energy-integrated
		 * fission source of mesh cells */
		eps = 0;
		for (int i = 0; i < cw * ch; i++)
		{
		    meshCell = _mesh->getCells(i);
			new_flux[i] = 0;
			/* integrates over energy */
			for (int e = 0; e < ng; e++)
			{
				xs = meshCell->getNuSigmaF()[e];
				new_flux[i] += xs * meshCell->getNewFlux()[e];
			} 
			eps += pow(new_flux[i] / old_flux[i] - 1.0, 2);
		}
		eps = pow(eps, 0.5);
	  
		for (int i = 0; i < cw * ch; i++)
			old_flux[i] = new_flux[i];
		

		log_printf(NORMAL, " %d-th LOO iteration k = %f, eps = %e", 
				   iter, _keff, eps);
		log_printf(ACTIVE, "  fission source = %f, abs source = %f", 
				   fis_tot, abs_tot);

		/* If LOO iterative solver converges */
		/* FIXME: previous version forces LOO to run at least 5 times */
		if (eps < _l2_norm_conv_thresh)
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
	
	log_printf(WARNING, "Keff not converging after %d iter", 
			   max_outer);

	/* Cleaning up; FIXME: more cleaning */
	for (int i = 0; i < cw * ch; i++)
	{
		delete[] sum_quad_flux[i];
		delete[] quad_xs[i];
		delete[] ratio[i];
		delete[] expo[i];
		delete[] tau[i];
		delete[] new_src[i];
		delete[] new_quad_src[i];
	}
	delete[] sum_quad_flux;
	delete[] quad_xs;
	delete[] ratio;
	delete[] tau;
	delete[] new_src;
	delete[] new_quad_src;


	return _keff;
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @param A matrix
 * @param M matrix
 * @param old flux vector
 * @param solve methed - either DIFFUSION or CMFD
 * @return petsc error indicator
 */
int Cmfd::constructAMPhi(Mat A, Mat M, Vec phi_old, solveType solveMethod){

	/* initialized variables */
	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int ng = NUM_ENERGY_GROUPS;
	MeshCell* meshCell;
	int petsc_err = 0;
	PetscInt indice1, indice2;
	PetscScalar value;

	/* if single group, set ng (number of groups) to 1 */
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
				indice1 = (PetscInt) ((y * cw + x) * ng + e);
				value = (PetscScalar) meshCell->getOldFlux()[e];
				petsc_err = VecSetValues(phi_old, 1, &indice1, 
										 &value, INSERT_VALUES);
				CHKERRQ(petsc_err);

				/* diagonal - A */

				/* add absorption term to diagonal */
				value = (PetscScalar) 
					meshCell->getSigmaA()[e] * meshCell->getVolume();
				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* add outscattering term to diagonal */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = meshCell->getSigmaS()[e*ng + g] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng+e;
						petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
						CHKERRQ(petsc_err);
					}
				}

				/* add fission terms to diagonal in M */
				for (int g = 0; g < ng; g++){
					value = meshCell->getChi()[e] * meshCell->getNuSigmaF()[g] * meshCell->getVolume();
					indice1 = (y*cw + x)*ng+e;
					indice2 = (y*cw + x)*ng + g;
					petsc_err = MatSetValues(M, 1, &indice1, 1, &indice2, &value, INSERT_VALUES);
				CHKERRQ(petsc_err);
				}


				/* add in scattering terms to off diagonals in A */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = - meshCell->getSigmaS()[g*ng + e] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng + g;
						petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
						CHKERRQ(petsc_err);
					}
				}

				/* RIGHT SURFACE */

				/* set transport term on diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(2)->getDHat()[e]      - meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(2)->getDDif()[e] * meshCell->getHeight();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1,1 , &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* set transport terms on off diagonals */
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

				/* LEFT SURFACE */

				/* set transport term on diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(0)->getDHat()[e]      + meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(0)->getDDif()[e] * meshCell->getHeight();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* set transport terms on off diagonals */
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

				/* BOTTOM SURFACE */

				/* set transport term on diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(1)->getDHat()[e]      - meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(1)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* set transport terms on off diagonals */
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

				/* TOP SURFACE */

				/* set transport term on diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(3)->getDHat()[e]      + meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 =  (y*cw + x)*ng + e;
				petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				CHKERRQ(petsc_err);

				/* set transport terms on off diagonals */
				if (y != 0){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(3)->getDHat()[e] - meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

					indice1 = (y*cw + x)*ng + e;
					indice2 = ((y-1)*cw + x)*ng + e;
					petsc_err = MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
					CHKERRQ(petsc_err);
				}
			}
		}
	}

	return petsc_err;
}


/* Update the MOC flux in each FSR
 * @param MOC iteration number
 */
void Cmfd::updateMOCFlux(int iteration){

	log_printf(INFO, "Updating MOC flux...");

	/* initialize variables */
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	double old_flux, new_flux;
	double* flux;
	int cw = _mesh->getCellWidth();
	int ch = _mesh->getCellHeight();
	int ng = NUM_ENERGY_GROUPS;
	int i, e;
	std::vector<int>::iterator iter;


	/* loop over mesh cells */
#if USE_OPENMP
#pragma omp parallel for private(meshCell, i, e, \
    new_flux, old_flux, flux, fsr, iter)
#endif 
	for (i = 0; i < cw * ch; i++){

		/* get pointer to current mesh cell */
		meshCell = _mesh->getCells(i);

		/* loop over groups */
		for (e = 0; e < ng; e++){

			/* get the old and new meshCell flux */
			if (_mesh->getMultigroup() == true){
				old_flux = meshCell->getOldFlux()[e];
				new_flux = meshCell->getNewFlux()[e];
			}
			else {
				old_flux = meshCell->getOldFlux()[0];
				new_flux = meshCell->getNewFlux()[0];
			}

			if (e == NUM_ENERGY_GROUPS - 1)
			{
				log_printf(INFO, "Update flux in Cell: %i,"
						   " old =  %e, new = %e, new/old = %f", 
						   i, old_flux, new_flux, new_flux / old_flux);
			}

			/* loop over FRSs in mesh cell */
			for (iter = meshCell->getFSRs()->begin(); 
				 iter != meshCell->getFSRs()->end(); ++iter) {
				fsr = &_flat_source_regions[*iter];

				/* get fsr flux */
				flux = fsr->getFlux();


				/* set new flux in FSR */
				flux[e] = new_flux / old_flux * flux[e];
			}
		}
	}
}



double Cmfd::computeDiffCorrect(double d, double h){

	if (_use_diffusion_correction == false)
	{
		return 1.0;
	}
	else
	{
		/* compute correction - F */
		double alpha, mu, expon;
		double rho, F;
		rho = 0.0;
		for (int p = 0; p < NUM_POLAR_ANGLES; p++){
			mu = std::cos(std::asin(_quad->getSinTheta(p)));
			expon = exp(- h / (3 * d * mu));
			alpha = (1 + expon) / (1 - expon) - 2 * mu / h;
			rho += mu * _quad->getWeight(p) * alpha;
		}

		F = 1 + h * rho / (2 * d);
		return F;
	}

}

/* compute the L2 norm of consecutive fission sources
 * @retun L2 norm
 */
int Cmfd::fisSourceNorm(Vec snew, int iter){

	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int ng = NUM_ENERGY_GROUPS;
	int petsc_err = 0;
	PetscScalar *old_source;
	PetscScalar *new_source;
	petsc_err = VecGetArray(_source_old, &old_source);
	petsc_err = VecGetArray(snew, &new_source);
	CHKERRQ(petsc_err);

	_l2_norm = 0.0;
	for (int i = 0; i < cw*ch; i++)
	{
		for (int e = 0; e < ng; e++)
		{
			if (new_source[i*ng+e] != 0.0)
				_l2_norm += pow(new_source[i * ng + e] / old_source[i * ng + e]
								- 1.0, 2);
		}
	}
	_l2_norm = pow(_l2_norm, 0.5);

	petsc_err = VecRestoreArray(_source_old, &old_source);
	petsc_err = VecRestoreArray(snew, &new_source);
	CHKERRQ(petsc_err);

	std::ofstream logfile;
	std::stringstream string;
	string << "l2_norm_" << (_num_azim*2) << "_" <<  _spacing << ".txt";
	std::string title_str = string.str();

	/* Write the message to the output file */
	if (iter == 1){
		logfile.open(title_str.c_str(), std::fstream::trunc);
		logfile << "iteration, l2_norm, keff" << std::endl;
	}
	else{
		logfile.open(title_str.c_str(), std::ios::app);
		logfile << iter << " " << _l2_norm << " " << _keff << std::endl;
	}

    logfile.close();

	return petsc_err;
} 

double Cmfd::getL2Norm(){
	return _l2_norm;
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

