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
 * Solver constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the trackgenerator
 */
Cmfd::Cmfd(Geometry* geom, Plotter* plotter, Mesh* mesh, bool updateFlux) {
	_geom = geom;
	_plotter = plotter;
	_mesh = mesh;
	_update_flux = updateFlux;
	_quad = new Quadrature(TABUCHI);

}

/**
 * cmfd Destructor clears all memory
 */
Cmfd::~Cmfd() {
}


/* compute the cross section for all MeshCells in the Mesh */
void Cmfd::computeXS(FlatSourceRegion* fsrs){

	_mesh->splitCorners();
	_flat_source_regions = fsrs;

	/* initialize variables */
	double volume, flux, abs, tot, nu_fis, chi;
	double* scat;
	double abs_tally_cell, nu_fis_tally_cell, dif_tally_cell, rxn_tally_cell, vol_tally_cell, tot_tally_cell;
	double nu_fis_tally, dif_tally, rxn_tally, abs_tally, tot_tally;
	double scat_tally_cell[NUM_ENERGY_GROUPS];

	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	Material* material;

	/* loop over mesh cells */
	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++){
		meshCell = _mesh->getCells(i);

		if (_mesh->getMultigroup() == false){
			abs_tally = 0.0;
			nu_fis_tally = 0.0;
			dif_tally = 0.0;
			tot_tally = 0.0;
			rxn_tally = 0.0;
		}

		/* loop over energy groups */
		std::vector<int>::iterator iter;
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			abs_tally_cell = 0;
			nu_fis_tally_cell = 0;
			dif_tally_cell = 0;
			rxn_tally_cell = 0;
			vol_tally_cell = 0;
			tot_tally_cell = 0;

			for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
				scat_tally_cell[g] = 0;
			}

			/* loop over FSRs in mesh cell */
			for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter){
				fsr = &_flat_source_regions[*iter];

				/* Gets FSR specific data. */
				material = fsr->getMaterial();
				chi = material->getChi()[e];
				volume = fsr->getVolume();
				flux = fsr->getFlux()[e];
				abs = material->getSigmaA()[e];
				tot = material->getSigmaT()[e];
				nu_fis = material->getNuSigmaF()[e];
				scat = material->getSigmaS();

				abs_tally_cell += abs * flux * volume;
				tot_tally_cell += tot * flux * volume;
				dif_tally_cell += flux  * volume / (3.0 * tot);
				nu_fis_tally_cell += nu_fis * flux * volume;
				rxn_tally_cell += flux * volume;
				vol_tally_cell += volume;


				for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
					scat_tally_cell[g] += scat[g*NUM_ENERGY_GROUPS + e] * flux * volume;
				}


				if (chi >= meshCell->getChi()[e]){
					meshCell->setChi(chi,e);
				}
			}

			if (_mesh->getMultigroup() == true){
				meshCell->setVolume(vol_tally_cell);
				meshCell->setSigmaA(abs_tally_cell / rxn_tally_cell, e);
				meshCell->setSigmaT(tot_tally_cell / rxn_tally_cell, e);
				meshCell->setNuSigmaF(nu_fis_tally_cell / rxn_tally_cell, e);
				meshCell->setDiffusivity(dif_tally_cell / rxn_tally_cell, e);
				meshCell->setOldFlux(rxn_tally_cell / vol_tally_cell, e);

				for (int g = 0; g < NUM_ENERGY_GROUPS; g++){
					meshCell->setSigmaS(scat_tally_cell[g] / rxn_tally_cell,e,g);
				}
			}
			else{
				abs_tally += abs_tally_cell;
				tot_tally += tot_tally_cell;
				nu_fis_tally += nu_fis_tally_cell;
				dif_tally += dif_tally_cell;
				rxn_tally += rxn_tally_cell;
			}
		}


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

	/* compute fis and abs rates based on tallied fluxes an XSs */
	double fis_tot = 0;
	double abs_tot = 0;
	for (int i = 0; i < _mesh->getCellWidth() * _mesh->getCellHeight(); i++){
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++){
			fis_tot += meshCell->getNuSigmaF()[e]*meshCell->getOldFlux()[e]*meshCell->getVolume();
			abs_tot += meshCell->getSigmaA()[e]*meshCell->getOldFlux()[e]*meshCell->getVolume();
		}
	}

	/* print keff based on nu_fis / abs */
	log_printf(DEBUG, "fission rate / abs rate: %f", fis_tot / abs_tot);
}


/* compute the xs for all MeshCells in the Mesh */
void Cmfd::computeDs(){

	/* initialize variables */
	double d, d_next, d_hat, d_tilde, current, flux, flux_next, f = 1, f_next = 1;
	MeshCell* meshCell;
	MeshCell* meshCellNext;
	int ng = 1;

	if (_mesh->getMultigroup() == true){
		ng = NUM_ENERGY_GROUPS;
	}

	/* set cell width and height */
	int cell_height = _mesh->getCellHeight();
	int cell_width = _mesh->getCellWidth();

	/* loop over mesh cells in y direction */
	for (int y = 0; y < cell_height; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cell_width; x++){

			/* get mesh cell */
			meshCell = _mesh->getCells(y*cell_width + x);

			for (int e = 0; e < ng; e++){

				/* get diffusivity and flux for mesh cell */
				d = meshCell->getDiffusivity()[e];
				flux = meshCell->getOldFlux()[e];

				/* get diffusion correction term for meshCell */
				f = computeDiffCorrect(d, meshCell->getWidth());

				/* LEFT */
				/* if cell on left side, set d_hat and d_tilde to 0 */
				if (x == 0){
					if (_mesh->getBoundary(0) == REFLECTIVE){
						d_hat = 0.0;
						d_tilde = 0.0;
						current = 0.0;

						/* set d_dif */
						meshCell->getMeshSurfaces(0)->setDDif(0.0, e);
					}
					else if (_mesh->getBoundary(0) == VACUUM){
						current = - meshCell->getMeshSurfaces(0)->getCurrent(e);

						/* set d_dif */
						meshCell->getMeshSurfaces(0)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);
						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
						d_tilde = - (d_hat * flux + current / meshCell->getHeight()) / flux;
					}
				}
				else{
					/* get mesh cell to left */
					meshCellNext = _mesh->getCells(y*cell_width + x - 1);
					d_next = meshCellNext->getDiffusivity()[e];
					flux_next = meshCellNext->getOldFlux()[e];

					/* set d_dif */
					meshCell->getMeshSurfaces(0)->setDDif(2.0 * d * d_next / (meshCell->getWidth() * d + meshCellNext->getWidth() * d_next), e);

					/* get diffusion correction term for meshCellNext */
					f_next = computeDiffCorrect(d_next, meshCellNext->getWidth());

					/* compute d_hat */
					d_hat = 2.0 * d*f * d_next*f_next / (meshCell->getWidth() * d*f + meshCellNext->getWidth() * d_next*f_next);

					/* get net outward current across surface */
					current = 0.0;

					/* increment current by outward current on next cell's right side */
					current += meshCellNext->getMeshSurfaces(2)->getCurrent(e);

					/* decrement current by outward current on left side */
					current -= meshCell->getMeshSurfaces(0)->getCurrent(e);

					/* compute d_tilde */
					d_tilde = -(d_hat * (flux - flux_next) + current  / meshCell->getHeight()) / (flux_next + flux);

				}

				/* if abs(d_tilde) > abs(d_hat) -> make them equal in magnitude */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(DEBUG, "correcting Ds");

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux*meshCell->getHeight());
						d_tilde = - current/(2*flux*meshCell->getHeight());
					}
					else{
						d_hat   = current/(2*flux_next*meshCell->getHeight());
						d_tilde = - current/(2*flux_next*meshCell->getHeight());
					}
				}

				log_printf(DEBUG, "cell: %i, group: %i, side:   LEFT, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

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

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(DEBUG, "correcting Ds");

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

				log_printf(DEBUG, "cell: %i, group: %i, side: BOTTOM, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

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
						meshCell->getMeshSurfaces(2)->setDDif(2 * d / meshCell->getHeight() / (1 + 4 * d / meshCell->getHeight()), e);

						d_hat = 2 * d*f / meshCell->getHeight() / (1 + 4 * d*f / meshCell->getHeight());
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

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(DEBUG, "correcting Ds");

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

				log_printf(DEBUG, "cell: %i, group: %i, side:  RIGHT, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);


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

				/* if abs(d_tilde) > abs(d_hat) -> make them equal to each other */
				if (fabs(d_tilde) > fabs(d_hat)){
					log_printf(DEBUG, "correcting Ds");

					/* d_tilde is positive */
					if (1 - fabs(d_tilde)/d_tilde < 1e-8){
						d_hat   = - current/(2*flux*meshCell->getHeight());
						d_tilde = - current/(2*flux*meshCell->getHeight());
					}
					else{
						d_hat   = current/(2*flux_next*meshCell->getHeight());
						d_tilde = - current/(2*flux_next*meshCell->getHeight());
					}
				}

				log_printf(DEBUG, "cell: %i, group: %i, side:    TOP, current: %f, dhat: %f, dtilde: %f", y*cell_width + x, e, current, d_hat, d_tilde);

				/* set d_hat and d_tilde */
				meshCell->getMeshSurfaces(3)->setDHat(d_hat, e);
				meshCell->getMeshSurfaces(3)->setDTilde(d_tilde, e);
			}
		}
	}

}


/*
 * compute the flux in each mesh cell using power iteration with Petsc's GMRES numerical inversion
 */
double Cmfd::computeCMFDFluxPower(solveType solveMethod, int moc_iter){

	log_printf(NORMAL, "Running diffusion solver...");

	MeshCell* meshCell;
	double keff;
	int ng = NUM_ENERGY_GROUPS;
	if (_mesh->getMultigroup() == false){
		ng = 1;
	}

	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int petsc_err;

	Mat A, M;
	Vec phi_old;
	PetscInt size1, size2;

	size1 = ch*cw*ng;
	size2 = 4 + ng;
	petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &A);
	size2 = ng;
	petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &M);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, size1, &phi_old);
	CHKERRQ(petsc_err);

	/* construct A matrix, M matrix, and flux vector */
	constructAMPhi(A, M, phi_old, solveMethod);

	int max_outer = 1000;
	PetscScalar sumold, sumnew, scale_val, eps;
	double criteria = 1e-8;

	Vec sold, snew, phi_new, res;
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &sold);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &snew);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &phi_new);
	petsc_err = VecCreateSeq(PETSC_COMM_WORLD, ch*cw*ng, &res);
	CHKERRQ(petsc_err);

	size1 = 1;
	petsc_err = VecSet(phi_new, size1);
	CHKERRQ(petsc_err);

	/* assembly vectors and matrices */
	petsc_err = VecAssemblyBegin(phi_old);
	petsc_err = VecAssemblyEnd(phi_old);
	petsc_err = VecAssemblyBegin(phi_new);
	petsc_err = VecAssemblyEnd(phi_new);
	petsc_err = VecAssemblyBegin(sold);
	petsc_err = VecAssemblyEnd(sold);
	petsc_err = VecAssemblyBegin(snew);
	petsc_err = VecAssemblyEnd(snew);
	petsc_err = VecAssemblyBegin(res);
	petsc_err = VecAssemblyEnd(res);
	CHKERRQ(petsc_err);

	petsc_err = VecCopy(phi_old, phi_new);

	petsc_err = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
	petsc_err = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
	CHKERRQ(petsc_err);

	/* initialize KSP solver */
	KSP ksp;
	petsc_err = KSPCreate(PETSC_COMM_WORLD, &ksp);
	petsc_err = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	petsc_err = KSPSetType(ksp, KSPGMRES);
	petsc_err = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	petsc_err = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);
	petsc_err = KSPSetUp(ksp);
	petsc_err = KSPSetFromOptions(ksp);
	CHKERRQ(petsc_err);

	petsc_err = MatMult(M, phi_new, sold);
	petsc_err = VecSum(sold, &sumold);
	scale_val = (cw * ch * ng) / sumold;
	petsc_err = VecScale(sold, scale_val);
	CHKERRQ(petsc_err);

	sumold = cw * ch * ng;
	int iter = 0;

	for (iter = 0; iter < max_outer; iter++){

		petsc_err = KSPSolve(ksp, sold, phi_new);
		petsc_err = MatMult(M, phi_new, snew);
		petsc_err = VecSum(snew, &sumnew);
		CHKERRQ(petsc_err);

		keff = sumnew / sumold;
		petsc_err = VecScale(sold, keff);
		scale_val = 1e-15;
		petsc_err = VecShift(snew, scale_val);
		petsc_err = VecPointwiseDivide(res, sold, snew);
		scale_val = -1;
		petsc_err = VecShift(res, scale_val);
		CHKERRQ(petsc_err);

		petsc_err = VecNorm(res, NORM_2, &eps);

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

	log_printf(NORMAL, "Diffusion solver iter: %i, keff: %f", iter, keff);

	petsc_err = VecSum(phi_new, &sumnew);
	scale_val = sumnew / (cw*ch*ng);
	petsc_err = VecScale(phi_new, scale_val);
	petsc_err = VecSum(phi_old, &sumold);
	scale_val = sumold / (cw*ch*ng);
	petsc_err = VecScale(phi_old, scale_val);
	CHKERRQ(petsc_err);

	PetscScalar *old_phi;
	PetscScalar *new_phi;

	petsc_err = VecGetArray(phi_old, &old_phi);
	petsc_err = VecGetArray(phi_new, &new_phi);
	CHKERRQ(petsc_err);

	for (int i = 0; i < cw*ch; i++){
		meshCell = _mesh->getCells(i);
		for (int e = 0; e < ng; e++){
			meshCell->setOldFlux(double(old_phi[i*ng + e]), e);
			meshCell->setNewFlux(double(new_phi[i*ng + e]), e);
		}
	}

	petsc_err = VecRestoreArray(phi_old, &old_phi);
	petsc_err = VecRestoreArray(phi_new, &new_phi);
	CHKERRQ(petsc_err);

	/* destroy matrices and vectors */
	petsc_err = MatDestroy(&A);
	petsc_err = MatDestroy(&M);
	petsc_err = VecDestroy(&phi_old);
	petsc_err = VecDestroy(&phi_new);
	petsc_err = VecDestroy(&snew);
	petsc_err = VecDestroy(&sold);
	petsc_err = VecDestroy(&res);
	CHKERRQ(petsc_err);

	if (_update_flux == true && solveMethod == CMFD){
		updateMOCFlux();
	}

	std::string string;

	if (solveMethod == DIFFUSION){
		if (_plotter->plotDiffusion() == true){
			string = "diff";
			_plotter->plotCMFDflux(_mesh, string, moc_iter);
		}
	}

	/* compute the total surface currents */
	if (_mesh->getMultigroup() == false){
		_mesh->computeTotCurrents();
	}

	if (_plotter->plotKeff()){
		_plotter->plotCMFDKeff(_mesh, moc_iter);
	}

	if (_plotter->plotCurrent()){
		string = "cmfd";
		_plotter->plotCMFDflux(_mesh, string, moc_iter);
	}

	if (solveMethod == CMFD){
		_mesh->setKeffCMFD(keff, moc_iter);
	}

	return keff;
}


void Cmfd::constructAMPhi(Mat A, Mat M, Vec phi_old, solveType solveMethod){

	int ch = _mesh->getCellHeight();
	int cw = _mesh->getCellWidth();
	int ng;
	MeshCell* meshCell;

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

				VecSetValues(phi_old, 1, &indice1, &meshCell->getOldFlux()[e], INSERT_VALUES);

				/* diagonal - A */
				value = meshCell->getSigmaA()[e] * meshCell->getVolume();
				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);

				/* scattering out */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = meshCell->getSigmaS()[e*ng + g] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng+e;
						MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
					}
				}


				/* diagonal - M */
				for (int g = 0; g < ng; g++){
					value = meshCell->getChi()[e] * meshCell->getNuSigmaF()[g] * meshCell->getVolume();
					indice1 = (y*cw + x)*ng+e;
					indice2 = (y*cw + x)*ng + g;
					MatSetValues(M, 1, &indice1, 1, &indice2, &value, INSERT_VALUES);
				}


				/* scattering in */
				for (int g = 0; g < ng; g++){
					if (e != g){
						value = - meshCell->getSigmaS()[g*ng + e] * meshCell->getVolume();
						indice1 = (y*cw + x)*ng+e;
						indice2 = (y*cw + x)*ng + g;
						MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
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
				MatSetValues(A, 1, &indice1,1 , &indice2, &value, ADD_VALUES);

				/* transport in */
				if (x != cw - 1){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(2)->getDHat()[e] + meshCell->getMeshSurfaces(2)->getDTilde()[e]) * meshCell->getHeight();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(2)->getDDif()[e] * meshCell->getHeight();

					indice1 = (y*cw + x)*ng + e;
					indice2 = (y*cw + x + 1)*ng + e;
					MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				}

				/* LEFT */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(0)->getDHat()[e]      + meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(0)->getDDif()[e] * meshCell->getHeight();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);

				/* transport in */
				if (x != 0){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(0)->getDHat()[e] - meshCell->getMeshSurfaces(0)->getDTilde()[e]) * meshCell->getHeight();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(0)->getDDif()[e] * meshCell->getHeight();

					indice1 = (y*cw + x)*ng + e;
					indice2 = (y*cw + x - 1)*ng + e;
					MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				}

				/* BELOW */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(1)->getDHat()[e]      - meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(1)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 = (y*cw + x)*ng + e;
				MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);

				/* transport in */
				if (y != ch - 1){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(1)->getDHat()[e] + meshCell->getMeshSurfaces(1)->getDTilde()[e]) * meshCell->getWidth();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(1)->getDDif()[e] * meshCell->getWidth();

					indice1 = (y*cw + x)*ng + e;
					indice2 = ((y+1)*cw + x)*ng + e;
					MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				}

				/* ABOVE */
				/* diagonal */
				if (solveMethod == CMFD)
					value = (meshCell->getMeshSurfaces(3)->getDHat()[e]      + meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
				else if (solveMethod == DIFFUSION)
					value = meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

				indice1 = (y*cw + x)*ng + e;
				indice2 =  (y*cw + x)*ng + e;
				MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);

				/* transport in */
				if (y != 0){
					if (solveMethod == CMFD)
						value = - (meshCell->getMeshSurfaces(3)->getDHat()[e] - meshCell->getMeshSurfaces(3)->getDTilde()[e]) * meshCell->getWidth();
					else if (solveMethod == DIFFUSION)
						value = - meshCell->getMeshSurfaces(3)->getDDif()[e] * meshCell->getWidth();

					indice1 = (y*cw + x)*ng + e;
					indice2 = ((y-1)*cw + x)*ng + e;
					MatSetValues(A, 1, &indice1, 1, &indice2, &value, ADD_VALUES);
				}
			}
		}
	}
}





/* update the MOC flux in each FSR and track fluxes at the boundaries */
void Cmfd::updateMOCFlux(){

	log_printf(NORMAL, "Updating MOC flux...");

	/* initialize variables */
	MeshCell* meshCell;
	FlatSourceRegion* fsr;
	double old_flux, new_flux, fsr_new_flux;
	double* flux;

	int cw = _mesh->getCellWidth();
	int ch = _mesh->getCellHeight();
	int ng = NUM_ENERGY_GROUPS;

	/* loop over mesh cells */
	for (int i = 0; i < cw * ch; i++){

		/* get mesh cell and flux values */
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

			log_printf(DEBUG, "Updating flux in meshCell: %i, flux ratio: %f", i, new_flux / old_flux);

			/* loop over FRSs in mesh cell */
			std::vector<int>::iterator iter;
			for (iter = meshCell->getFSRs()->begin(); iter != meshCell->getFSRs()->end(); ++iter) {
				fsr = &_flat_source_regions[*iter];
				/* get fsr flux */
				flux = fsr->getFlux();
				fsr_new_flux = new_flux / old_flux * flux[e];

				/* set new flux in FSR */
				fsr->setFlux(e, fsr_new_flux);
			}
		}
	}

//	/* update track fluxes */
//	Track* track;
//	int num_segments;
//	std::vector<segment*> segments;
//	segment* segment;
//	double* polar_fluxes;
//	int e, p, pe;
//	double current, d_hat, d_tilde, cell_flux;
//	int surfId;
//
//
//	for (int j = 0; j < _num_azim; j++){
//
//		/* Loop over all tracks for this azimuthal angles */
//		for (int k = 0; k < _num_tracks[j]; k++) {
//
//			/* Initialize local pointers to important data structures */
//			track = &_tracks[j][k];
//			polar_fluxes = track->getPolarFluxes();
//
//			/* get first segment */
//			segments = track->getSegments();
//			num_segments = track->getNumSegments();
//			segment = segments.at(0);
//
//			/* set polar angle * energy group to 0 */
//			pe = 0;
//
//			/* loop over energy groups */
//			for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//				/* loop over polar angles */
//				for (p = 0; p < NUM_POLAR_ANGLES; p++){
//
//					surfId = segment->_mesh_surface_bwd->getId();
//
//					if (surfId == 0 || surfId == 1 || surfId == 2 || surfId == 3){
//						current = segment->_mesh_surface_bwd->getCurrent(e);
//						d_hat = segment->_mesh_surface_bwd->getDHat()[e];
//						d_tilde = segment->_mesh_surface_bwd->getDTilde()[e];
//						meshCell = mesh->getCells(segment->_mesh_surface_bwd->getCellId());
//						cell_flux = meshCell->getNewFlux()[e];
//
//						polar_fluxes[pe] = fabs(polar_fluxes[pe] / current * (d_hat - d_tilde) * cell_flux);
//
//						pe++;
//					}
//				}
//			}
//
//			/* get last segment */
//			segment = segments.at(num_segments - 1);
//
//
//			/* set polar angle * energy group to GRP_TIMES_ANG */
//			pe = GRP_TIMES_ANG;
//
//			/* loop over energy groups */
//			for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
//
//				/* loop over polar angles */
//				for (p = 0; p < NUM_POLAR_ANGLES; p++){
//
//					surfId = segment->_mesh_surface_fwd->getId();
//
//					if (surfId == 0 || surfId == 1 || surfId == 2 || surfId == 3){
//						current = segment->_mesh_surface_fwd->getCurrent(e);
//						d_hat = segment->_mesh_surface_fwd->getDHat()[e];
//						d_tilde = segment->_mesh_surface_fwd->getDTilde()[e];
//						meshCell = mesh->getCells(segment->_mesh_surface_fwd->getId());
//						cell_flux = meshCell->getNewFlux()[e];
//
//						polar_fluxes[pe] = fabs(polar_fluxes[pe] / current * (d_hat - d_tilde) * cell_flux);
//
//						pe++;
//					}
//				}
//			}
//		}
}



double Cmfd::computeDiffCorrect(double d, double h){

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





