/*
 * Transient.cpp
 *
 *  Created on: Dec 15, 2012
 *      Author: samuelshaner
 */

#include "Transient.h"

/**
 * Transient constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the trackgenerator
 */
Transient::Transient(Geometry* geom, Cmfd* cmfd, Mesh* mesh, Solver* solver, Plotter* plotter, double tMax, double dtOuter, double dtInner){

	_cmfd = cmfd;
	_geom = geom;
	_solver = solver;
	_plotter = plotter;

	/* copy geometry mesh to transient mesh */
	_mesh = new Mesh;

	_lambda[0] = 0.0654;
	_lambda[1] = 1.35;

	_amp[0] = 0.0;
	_amp[1] = 0.0;

	/* set delayed neutron precursor yield for group 0 and 1 */
	_beta[0] = 0.0054;
	_beta[1] = 0.001087;

	/* set velocities of group 0 and 1 */
	_velocities[0] = 3.0e7;
	_velocities[1] = 3.0e5;

	_alpha = 3.83e-11;
	_gamma = 3.034e-3;
	_kappa = 3.204e-11;
	_nu = 2.43;
	_buckle = 1.0e-4;
	_temp_0 = 300;
	_power_init = 1e-6;
	_vol_core = 17550;

	_t_max = tMax;
	_dt_outer = dtOuter;
	_dt_inner = dtInner;
	_t_steps = floor(_t_max / _dt_outer);
	_mesh->setOldTime(0.0);
	_geom->getMesh()->setOldTime(0.0);

	/* make arrays to store temperatures and powers at each step */
	_temp = new double[_t_steps+1];
	_power = new double[_t_steps+1];
	_power_pke = new double[_t_steps+1];
	_power_ass = new double[_t_steps+1];
	_temp_max = new double[_t_steps+1];

	_sigma_A_1 = _geom->getMaterial(5)->getSigmaA()[1];
	_sigma_S_4 = _geom->getMaterial(5)->getSigmaS()[3];

	MeshCell* meshCell;
	for (int i = 0; i < _geom->getMesh()->getCellHeight()*_geom->getMesh()->getCellWidth(); i++){
		meshCell = _geom->getMesh()->getCells(i);
		meshCell->setTemp(300);

		/* give mesh cell a pointer to material */
		meshCell->setMaterial(_solver->getFSRs()[meshCell->getFSRStart()].getMaterial());
	}
}

/**
 * Transient Destructor clears all memory
 */
Transient::~Transient() {

}

arma::mat Transient::constructM(Mesh* mesh, arma::mat M){

	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int n_cells = cw*ch;

	/* loop over mesh cells in y direction */
	for (int y = 0; y < ch; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cw; x++){

			/* get mesh cell */
			meshCell = mesh->getCells(y*cw+x);

			/* fission - M */
			M(0, y*cw+x)           = meshCell->getNuSigmaF()[0] * meshCell->getVolume();
			M(1, n_cells + y*cw+x) = meshCell->getNuSigmaF()[1] * meshCell->getVolume();

		}
	}

	return M;
}


arma::mat Transient::constructN(Mesh* mesh, arma::mat N, arma::mat A, arma::mat M, arma::colvec shape){

	using namespace arma;

	double beta = _beta[0] + _beta[1];

	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int nc = cw*ch;
	double g1_mult = 0.0, g2_mult = 0.0;

//	for (int i = 0; i < nc; i++){
//		meshCell = mesh->getCells(i);
//		g1_mult += meshCell->getVolume() * shape(i);
//		g2_mult += meshCell->getVolume() * shape(nc + i);
//	}
//
//	g1_mult = _velocities[0] / g1_mult;
//	g2_mult = _velocities[1] / g2_mult;

	for (int i = 0; i < nc; i++){
		meshCell = mesh->getCells(i);
		g1_mult += (meshCell->getVolume() * shape(i)) / _velocities[0];
		g2_mult += (meshCell->getVolume() * shape(nc + i)) / _velocities[1];
	}

	g1_mult = 1.0 / g1_mult;
	g2_mult = 1.0 / g2_mult;


	colvec shape_1 = shape(span(0,nc-1));
	colvec shape_2 = shape(span(nc,2*nc-1));
	N(0,0) = sum(((1-beta)*M(0,span::all)/_keff_0 * shape  - sum(A(span(0,4),span::all) * shape))) * g1_mult;
	N(0,1) = sum(((1-beta)*M(1,span::all)/_keff_0 * shape)) * g1_mult;
	N(1,0) = sum(- A(5,span::all) * shape ) * g2_mult;
	N(1,1) = sum(- A(span(6,10),span::all) * shape ) * g2_mult;

	N(2,0) = sum(_beta[0] / _keff_0 * M(0,span::all) * shape);
	N(2,1) = sum(_beta[0] / _keff_0 * M(1,span::all) * shape);
	N(3,0) = sum(_beta[1] / _keff_0 * M(0,span::all) * shape);
	N(3,1) = sum(_beta[1] / _keff_0 * M(0,span::all) * shape);

	N(0,2) = _lambda[0]*g1_mult;
	N(0,3) = _lambda[1]*g1_mult;
	N(2,2) = - _lambda[0];
	N(3,3) = - _lambda[1];

	return N;
}

arma::colvec Transient::constructShape(arma::colvec shape, Mesh* mesh){

	MeshCell* meshCell;

	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int n_cells = cw*ch;

	for (int i = 0; i < n_cells; i++){
		meshCell = mesh->getCells(i);

		shape(i) = meshCell->getOldFlux()[0];
		shape(n_cells + i) = meshCell->getOldFlux()[1];
	}

	return shape;
}

arma::colvec Transient::constructB(arma::colvec B, arma::mat N, arma::mat M, arma::colvec shape){

	double fis = sum(M*shape);

//	B(0) = 1.0;
//	B(1) = 1.0;
//	B(2) = _beta[0] / (_lambda[0] * _keff_0) * fis;
//	B(3) = _beta[1] / (_lambda[1] * _keff_0) * fis;

	B(2) = _beta[0] / (_lambda[0] * _keff_0) * fis;
	B(3) = _beta[1] / (_lambda[1] * _keff_0) * fis;
	B(1) = (N(0,2)*B(2) + N(0,3)*B(3))/(N(0,0)*N(1,1)/N(1,0) - N(0,1));
	B(0) = -N(1,1)/N(1,0)*B(1);

	return B;
}


void Transient::solve(){

	log_printf(NORMAL, "Solving time dependent problem...");

	using namespace arma;

	/* initialize variables */
	int cw = _geom->getMesh()->getCellWidth();
	int ch = _geom->getMesh()->getCellHeight();
	double keff = 0.0;

	colvec shape = zeros(ch*cw*2);
	mat N = zeros(4,4);
	vec B = zeros(4);
	vec C = zeros(4);
	colvec resid = zeros(4);
	colvec D = zeros(4);
	mat A = zeros(11,ch*cw*2);
	mat M = zeros(2,ch*cw*2);

	/* get keff_0 */
	_keff_0 = _solver->computeKeff(3000);

	/* get initial matrices and vectors */
	A = constructA(_geom->getMesh(), A, CMFD);
	M = constructM(_geom->getMesh(), M);
	shape = constructShape(shape, _geom->getMesh());
	N = constructN(_geom->getMesh(), N, A, M, shape);
	B = constructB(B, N, M, shape);
	C = B;

	_temp[0] = 300;
	_power[0] = 1e-6;
	_amp[0] = B(0);
	_amp[1] = B(1);
	_precursors[0] = B(2);
	_precursors[1] = B(3);
	_power_pke[0] = computePower(_geom->getMesh(), B);
	_power_factor = _power_init / _power_pke[0];
	log_printf(NORMAL, "power conversion factor: %f", _power_factor);
	_power[0] = _power_pke[0] * _power_factor;

	int iteration = 0;
	for (double time = 0.0; time < _t_max; time += _dt_outer){
		iteration++;

		/* move blade */
//		moveBlade(_geom->getMesh(), time);

		/* run MOC solver */
		keff = _solver->computeKeff(3000);

		/* get new shape, A, and N */
		N.zeros();
		shape.zeros();
		A.zeros();
		A = constructA(_geom->getMesh(), A, CMFD);
		shape = constructShape(shape, _geom->getMesh());
		N = constructN(_geom->getMesh(), N, A, M, shape);
		C = B;

		std::cout << N << endl;
		std::cout << B << endl;

		/* step forward in time - BWD Euler */
		int counter = 0;
		for (double t = time; t < time + _dt_outer; t += _dt_inner){
			for (int inner = 0; inner < 25; inner++){
				D = B + _dt_inner * N * C;
				resid = D - C;
				C = D;
				if (fabs(sum(resid)) < 1e-8){
					break;
				}
			}
			B = C;

			counter++;
			if (counter == 1e2){
//				moveBlade(_geom->getMesh(), t);
//				tempFeedback(_geom->getMesh(), B, t);
				A.zeros();
				N.zeros();
				A = constructA(_geom->getMesh(), A, CMFD);
				N = constructN(_geom->getMesh(), N, A, M, shape);
				counter = 0;
			}
		}

		/* print new B */
		std::cout << B << endl;

		/* temperature feedback */
		tempFeedback(_geom->getMesh(), B, time + _dt_outer);
		_power[iteration]     = computePower(_geom->getMesh(), B) * _power_factor;
		_temp[iteration]      = computeTemp(_geom->getMesh(), B);
		_power_ass[iteration] = computePowerAss(_geom->getMesh(), B);
		_temp_max[iteration]  = computeTempMax(_geom->getMesh(), B);

		log_printf(NORMAL, "power ass: %f", _power_ass[iteration]);

		plotTemp(iteration, time+ _dt_outer);
		plotPower(iteration, time+ _dt_outer);
		plotTempMax(iteration, time+ _dt_outer);
		plotPowerAss(iteration, time+ _dt_outer);
		plotTempAss(_geom->getMesh());

		log_printf(NORMAL, "TIME: %f, POWER: %f, TEMP: %f", time, _power[iteration], _temp[iteration]);
	}

}



void Transient::plotTempAss(Mesh* mesh){

	log_printf(NORMAL, "plotting assembly temperatures...");

	int _bit_length_x = _plotter->getBitLengthX();
	int _bit_length_y = _plotter->getBitLengthY();

	/* set up bitMap */
	BitMap<double>* bitMap = new BitMap<double>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMap->color_type = SCALED;

	double x_global;
	double y_global;
	std::stringstream string;
	std::stringstream num;
	std::string title_str;

	/* find meshCell for each pixel */
	for (int y=0;y < _bit_length_y; y++){
		for (int x = 0; x < _bit_length_x; x++){
			x_global = _plotter->convertToGeometryX(x);
			y_global = _plotter->convertToGeometryY(y);
			bitMap->pixels[y * _bit_length_x + x] = mesh->getCells(mesh->findMeshCell(x_global, y_global))->getTemp();
		}
	}

	string.str("");
	string << "temp_ass";
	title_str = string.str();

	plot(bitMap, title_str, "png");

	deleteBitMap(bitMap);

}


/* change xs due to temperature changes and move control blade */
void Transient::tempFeedback(Mesh* mesh, arma::colvec B, double time_new){

	/* compute power multiplication factor */
	double cell_power, cell_temp;
	double abs_base, abs_new, scat_base;


	MeshCell* meshCell;
	FlatSourceRegion* fsrs;
	fsrs = _solver->getFSRs();
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int nc = cw*ch;

	/* loop over cells */
	for (int i = 0; i < nc; i++){
		meshCell = mesh->getCells(i);

		cell_power = 0.0;

		/* compute power in cell */
		cell_power += _alpha / _nu * meshCell->getNuSigmaF()[0]*meshCell->getOldFlux()[0]*meshCell->getVolume()*B(0);
		cell_power += _alpha / _nu * meshCell->getNuSigmaF()[1]*meshCell->getOldFlux()[1]*meshCell->getVolume()*B(1);

		/* get actual power in W/cc */
		cell_power = cell_power * _power_factor;

		/* get cell temp */
		cell_temp = meshCell->getTemp() + (time_new - mesh->getOldTime()) * cell_power;

		meshCell->setTemp(cell_temp);

		/* adjust abs xs for FSRs in that cell */
		for (int fsr_id = meshCell->getFSRStart(); fsr_id <= meshCell->getFSREnd(); fsr_id++){

			abs_base = fsrs[fsr_id].getMaterial()->getSigmaA()[0];
			scat_base = fsrs[fsr_id].getMaterial()->getSigmaS()[0];
			abs_new = abs_base * (1 + _gamma * (sqrt(cell_temp) - sqrt(300)));
			fsrs[fsr_id].setMatMultA(0, abs_new/abs_base);
			fsrs[fsr_id].setMatMult(0, (scat_base + abs_base - abs_new)/scat_base);
		}

		/* adjust abs and tot xs for meshcell */
		meshCell->setSigmaA(abs_new,0);
		meshCell->setSigmaS(scat_base + abs_base - abs_new,0,0);
	}

	mesh->setOldTime(time_new);

}


void Transient::moveBlade(Mesh* mesh, double time){

	double* sigmaA;
	double* sigmaS;
	sigmaA = _geom->getMaterial(5)->getSigmaA();
	sigmaS = _geom->getMaterial(5)->getSigmaS();

	/* move blade for global materials */
	if (time < 2.0){
		sigmaA[1] = _sigma_A_1 * (1 - (.0606184 * time));
		sigmaS[3] = _sigma_S_4 + _sigma_A_1 * (.0606184 * time);
	}
	else{
		sigmaA[1] = _sigma_A_1 * (1 - (.0606184 * 2));
		sigmaS[3] = _sigma_S_4 + _sigma_A_1 * (.0606184 * 2);
	}

	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int nc = cw*ch;

	/* move blade for local mesh */
	for (int i = 0; i < nc; i++){
		meshCell = mesh->getCells(i);

		if (meshCell->getMaterial()->getId() == 5){
			meshCell->setSigmaA(sigmaA[1], 1);
			meshCell->setSigmaS(sigmaS[3], 1,1);
		}
	}
}



/* copy MeshCell and MeshSurface values from mesh1 to mesh2 */
void Transient::copyMesh(Mesh* mesh1, Mesh* mesh2){

	MeshCell* meshCell1;
	MeshCell* meshCell2;
	int cw = mesh1->getCellWidth();
	int ch = mesh1->getCellHeight();

	for (int i = 0; i < cw*ch; i++){
		meshCell1 = mesh1->getCells(i);
		meshCell2 = mesh2->getCells(i);

		meshCell2->setVolume(meshCell1->getVolume());
		meshCell2->setFSRStart(meshCell1->getFSRStart());
		meshCell2->setFSREnd(meshCell1->getFSREnd());
		meshCell2->setTemp(meshCell1->getTemp());

		/* loop over energy groups */
		for (int e = 0; e < 2; e++){

			/* copy new sigA, sigT, sigS, nuSigF, diff, flux */
			meshCell2->setSigmaA(meshCell1->getSigmaA()[e],e);
			meshCell2->setSigmaT(meshCell1->getSigmaT()[e],e);
			meshCell2->setNuSigmaF(meshCell1->getNuSigmaF()[e],e);
			meshCell2->setDiffusivity(meshCell1->getDiffusivity()[e],e);
			meshCell2->setOldFlux(meshCell1->getOldFlux()[e],e);
			meshCell2->setNewFlux(meshCell1->getNewFlux()[e],e);

			for (int g = 0; g < 2; g++){
				meshCell2->setSigmaS(meshCell1->getSigmaS()[e*2+g],e,g);
			}

			/* loop over surfaces */
			for (int s = 0; s < 8; s++){

				/* copy dHat, dTilde, dDif */
				meshCell2->getMeshSurfaces(s)->setDHat(meshCell1->getMeshSurfaces(s)->getDHat()[e],e);
				meshCell2->getMeshSurfaces(s)->setDTilde(meshCell1->getMeshSurfaces(s)->getDTilde()[e],e);
				meshCell2->getMeshSurfaces(s)->setDDif(meshCell1->getMeshSurfaces(s)->getDDif()[e],e);

			}
		}
	}

}

/* compute core temperature */
double Transient::computeTemp(Mesh* mesh, arma::colvec B){

	double mesh_temp = 0.0;
	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int nc = cw*ch;

	/* loop over cells */
	for (int i = 0; i < nc; i++){
		meshCell = mesh->getCells(i);

		if (meshCell->getMaterial()->getId() < 6){
			mesh_temp += meshCell->getTemp() * meshCell->getVolume();
		}
	}

	mesh_temp = mesh_temp / _vol_core;

	return mesh_temp;
}


/* compute core temperature */
double Transient::computeTempMax(Mesh* mesh, arma::colvec B){

	double temp_max = 0.0;
	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int nc = cw*ch;

	/* loop over cells */
	for (int i = 0; i < nc; i++){
		meshCell = mesh->getCells(i);

		if (meshCell->getTemp() > temp_max){
			temp_max = meshCell->getTemp();
		}
	}

	return temp_max;
}



/* compute the power / cc */
double Transient::computePower(Mesh* mesh, arma::colvec B){

	MeshCell* meshCell;

	double power = 0.0;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();

	for (int i = 0; i < cw*ch; i++){
		meshCell = mesh->getCells(i);

		power += _kappa / _nu * meshCell->getNuSigmaF()[0]*meshCell->getOldFlux()[0]*meshCell->getVolume() * B(0);
		power += _kappa / _nu * meshCell->getNuSigmaF()[1]*meshCell->getOldFlux()[1]*meshCell->getVolume() * B(1);
	}

	power = power / _vol_core;

	return power;
}


/* compute the max normalized assembly power */
double Transient::computePowerAss(Mesh* mesh, arma::colvec B){

	MeshCell* meshCell;

	double power = 0.0, power_max = 0.0, power_cell = 0.0;;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();

	for (int i = 0; i < cw*ch; i++){
		meshCell = mesh->getCells(i);
		power_cell = 0.0;

		power_cell += _kappa / _nu * meshCell->getNuSigmaF()[0]*meshCell->getOldFlux()[0]*meshCell->getVolume() * B(0);
		power_cell += _kappa / _nu * meshCell->getNuSigmaF()[1]*meshCell->getOldFlux()[1]*meshCell->getVolume() * B(1);

		power += power_cell;

		if (power_cell < power_max){
			power_max = power_cell;
		}
	}

	/* max norm power */
	power_max = power_max / (power / (78*(225/meshCell->getVolume())));

	return power_max;
}



arma::mat Transient::constructA(Mesh* mesh, arma::mat A, solveType solveMethod){

	MeshCell* meshCell;
	int cw = mesh->getCellWidth();
	int ch = mesh->getCellHeight();
	int n_cells = cw*ch;

	/* loop over mesh cells in y direction */
	for (int y = 0; y < ch; y++){

		/* loop over mesh cells in x direction */
		for (int x = 0; x < cw; x++){

			/* get mesh cell */
			meshCell = mesh->getCells(y*cw+x);

			/* absorption - A */
			A(2,y*cw+x) += meshCell->getSigmaA()[0] * meshCell->getVolume();
			A(8, n_cells + y*cw+x) += meshCell->getSigmaA()[1] * meshCell->getVolume();

			/* down scattering - A */
			A(2,y*cw+x) += meshCell->getSigmaS()[1] * meshCell->getVolume();
			A(5, y*cw+x) -= meshCell->getSigmaS()[1] * meshCell->getVolume();


			/* RIGHT */
			/* diagonal */
			if (solveMethod == CMFD){
				A(2, y*cw+x)                     += (meshCell->getMeshSurfaces(2)->getDHat()[0]      - meshCell->getMeshSurfaces(2)->getDTilde()[0]) * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += (meshCell->getMeshSurfaces(2)->getDHat()[1]      - meshCell->getMeshSurfaces(2)->getDTilde()[1]) * meshCell->getHeight();			}
			else if (solveMethod == DIFFUSION){
				A(2, y*cw+x)                     += meshCell->getMeshSurfaces(2)->getDDif()[0] * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += meshCell->getMeshSurfaces(2)->getDDif()[1] * meshCell->getHeight();
			}

			/* transport in */
			if (x != cw - 1){
				if (solveMethod == CMFD){
					A(3, y*cw+x + 1)                     -= (meshCell->getMeshSurfaces(2)->getDHat()[0] + meshCell->getMeshSurfaces(2)->getDTilde()[0]) * meshCell->getHeight();
					A(9, n_cells + y*cw+x + 1) -= (meshCell->getMeshSurfaces(2)->getDHat()[1] + meshCell->getMeshSurfaces(2)->getDTilde()[1]) * meshCell->getHeight();
				}
				else if (solveMethod == DIFFUSION){
					A(3, y*cw+x + 1)                     -= meshCell->getMeshSurfaces(2)->getDDif()[0] * meshCell->getHeight();
					A(9, n_cells + y*cw+x + 1) -= meshCell->getMeshSurfaces(2)->getDDif()[1] * meshCell->getHeight();
				}
			}

			/* LEFT */
			/* diagonal */
			if (solveMethod == CMFD){
				A(2, y*cw+x)                     += (meshCell->getMeshSurfaces(0)->getDHat()[0]      - meshCell->getMeshSurfaces(0)->getDTilde()[0]) * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += (meshCell->getMeshSurfaces(0)->getDHat()[1]      - meshCell->getMeshSurfaces(0)->getDTilde()[1]) * meshCell->getHeight();			}
			else if (solveMethod == DIFFUSION){
				A(2, y*cw+x)                     += meshCell->getMeshSurfaces(0)->getDDif()[0] * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += meshCell->getMeshSurfaces(0)->getDDif()[1] * meshCell->getHeight();
			}

			if (x != 0){
				if (solveMethod == CMFD){
					A(1, y*cw+x - 1)                     -= (meshCell->getMeshSurfaces(0)->getDHat()[0] + meshCell->getMeshSurfaces(2)->getDTilde()[0]) * meshCell->getHeight();
					A(7, n_cells + y*cw+x - 1) -= (meshCell->getMeshSurfaces(0)->getDHat()[1] + meshCell->getMeshSurfaces(2)->getDTilde()[1]) * meshCell->getHeight();
				}
				else if (solveMethod == DIFFUSION){
					A(1, y*cw+x - 1)                     -= meshCell->getMeshSurfaces(0)->getDDif()[0] * meshCell->getHeight();
					A(7, n_cells + y*cw+x - 1) -= meshCell->getMeshSurfaces(0)->getDDif()[1] * meshCell->getHeight();
				}
			}

			/* BELOW */
			/* diagonal */
			if (solveMethod == CMFD){
				A(2, y*cw+x)                     += (meshCell->getMeshSurfaces(1)->getDHat()[0]      - meshCell->getMeshSurfaces(1)->getDTilde()[0]) * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += (meshCell->getMeshSurfaces(1)->getDHat()[1]      - meshCell->getMeshSurfaces(1)->getDTilde()[1]) * meshCell->getHeight();			}
			else if (solveMethod == DIFFUSION){
				A(2, y*cw+x)                     += meshCell->getMeshSurfaces(1)->getDDif()[0] * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += meshCell->getMeshSurfaces(1)->getDDif()[1] * meshCell->getHeight();
			}

			if (y != ch - 1){
				if (solveMethod == CMFD){
					A(4, y*cw+x + ch)                     -= (meshCell->getMeshSurfaces(1)->getDHat()[0] + meshCell->getMeshSurfaces(1)->getDTilde()[0]) * meshCell->getHeight();
					A(10, n_cells + y*cw+x + ch) -= (meshCell->getMeshSurfaces(1)->getDHat()[1] + meshCell->getMeshSurfaces(1)->getDTilde()[1]) * meshCell->getHeight();
				}
				else if (solveMethod == DIFFUSION){
					A(4, y*cw+x + ch)                     -= meshCell->getMeshSurfaces(1)->getDDif()[0] * meshCell->getHeight();
					A(10, n_cells + y*cw+x + ch) -= meshCell->getMeshSurfaces(1)->getDDif()[1] * meshCell->getHeight();
				}
			}

			/* ABOVE */
			/* diagonal */
			if (solveMethod == CMFD){
				A(2, y*cw+x)                     += (meshCell->getMeshSurfaces(3)->getDHat()[0]      - meshCell->getMeshSurfaces(3)->getDTilde()[0]) * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += (meshCell->getMeshSurfaces(3)->getDHat()[1]      - meshCell->getMeshSurfaces(3)->getDTilde()[1]) * meshCell->getHeight();			}
			else if (solveMethod == DIFFUSION){
				A(2, y*cw+x)                     += meshCell->getMeshSurfaces(3)->getDDif()[0] * meshCell->getHeight();
				A(8, n_cells + y*cw+x) += meshCell->getMeshSurfaces(3)->getDDif()[1] * meshCell->getHeight();
			}

			if (y != 0){
				if (solveMethod == CMFD){
					A(0, y*cw+x - ch)                     -= (meshCell->getMeshSurfaces(3)->getDHat()[0] + meshCell->getMeshSurfaces(3)->getDTilde()[0]) * meshCell->getHeight();
					A(6, n_cells + y*cw+x - ch) -= (meshCell->getMeshSurfaces(3)->getDHat()[1] + meshCell->getMeshSurfaces(3)->getDTilde()[1]) * meshCell->getHeight();
				}
				else if (solveMethod == DIFFUSION){
					A(0, y*cw+x - ch)                     -= meshCell->getMeshSurfaces(3)->getDDif()[0] * meshCell->getHeight();
					A(6, n_cells + y*cw+x - ch) -= meshCell->getMeshSurfaces(3)->getDDif()[1] * meshCell->getHeight();
				}
			}
		}
	}

	return A;
}

Mesh* Transient::getMesh(){
	return _mesh;
}


void Transient::plotTemp(int num_iter, double time){
	log_printf(NORMAL, "plotting Temp...");

	int _bit_length_x = _plotter->getBitLengthX();
	int _bit_length_y = _plotter->getBitLengthY();

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _bit_length_x;
	bitMap->geom_y = _bit_length_y;
	bitMap->color_type = BLACKWHITE;

	std::stringstream text_stream;
	std::string text;
	std::string text2;

	/* draw and label axes */
	drawLine(bitMap, _bit_length_x / 10, 10, _bit_length_x / 10, _bit_length_y - 10);
	drawLine(bitMap, 10, 9 * _bit_length_y / 10, _bit_length_x - 10, 9 * _bit_length_y / 10);
	text = "temp";
	drawText(bitMap, text, 30, _bit_length_y / 2 - 30);
	text = "time";
	drawText(bitMap, text, _bit_length_x / 2, _bit_length_y - 30);
	text_stream << time;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, 9.5 * _bit_length_x / 10, 9 * _bit_length_y / 10 + 20);


	/* create x axis scale */
	double temp_max = 0, temp_min = 2;
	for (int i = 1; i < num_iter; i++){
		if (_temp[i] > temp_max){
			temp_max = _temp[i];
		}
		if (_temp[i] > temp_max){
			temp_max = _temp[i];
		}
		if (_temp[i] < temp_min){
			temp_min = _temp[i];
		}
		if (_temp[i] < temp_min){
			temp_min = _temp[i];
		}
	}

	text_stream << temp_min;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, 9 * _bit_length_y / 10 - 10);
	text_stream << temp_max;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 10);
	double temp_avg = (temp_max - temp_min)/2 + temp_min;
	double temp_range = temp_max - temp_min;
	text_stream << temp_avg;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 2);

	text = "blue";
	text2 = "white";

	/* draw CMFD keff */
	for (int i = 1; i < num_iter; i++){
		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 2 - (_temp[i] - temp_avg) / (temp_range / 2) * (4 * _bit_length_y / 10), 3);
	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30, 3);
	text = "Temp";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 5);

	text = "temp";

	/* create filename with correct extension */
	plot(bitMap, text, "png");

	deleteBitMap(bitMap);

}

void Transient::plotPower(int num_iter, double time){
	log_printf(NORMAL, "plotting Power...");

	int _bit_length_x = _plotter->getBitLengthX();
	int _bit_length_y = _plotter->getBitLengthY();

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _bit_length_x;
	bitMap->geom_y = _bit_length_y;
	bitMap->color_type = BLACKWHITE;

	std::stringstream text_stream;
	std::string text;
	std::string text2;

	/* draw and label axes */
	drawLine(bitMap, _bit_length_x / 10, 10, _bit_length_x / 10, _bit_length_y - 10);
	drawLine(bitMap, 10, 9 * _bit_length_y / 10, _bit_length_x - 10, 9 * _bit_length_y / 10);
	text = "power";
	drawText(bitMap, text, 30, _bit_length_y / 2 - 30);
	text = "time";
	drawText(bitMap, text, _bit_length_x / 2, _bit_length_y - 30);
	text_stream << time;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, 9.5 * _bit_length_x / 10, 9 * _bit_length_y / 10 + 20);


	/* create x axis scale */
	double power_max = 0, power_min = 2;
	for (int i = 1; i < num_iter; i++){
		if (_power[i] > power_max){
			power_max = _power[i];
		}
		if (_power[i] > power_max){
			power_max = _power[i];
		}
		if (_power[i] < power_min){
			power_min = _power[i];
		}
		if (_power[i] < power_min){
			power_min = _power[i];
		}
	}

	text_stream << power_min;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, 9 * _bit_length_y / 10 - 10);
	text_stream << power_max;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 10);
	double power_avg = pow(10,(log10(power_max) + log10(power_min))/2);
	text_stream << power_avg;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 2);

	text = "blue";
	text2 = "white";

	/* draw CMFD keff */
	for (int i = 1; i < num_iter; i++){
		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 10 + (1 - log10(_power[i]/power_min) / log10(power_max / power_min)) * (8 * _bit_length_y / 10), 3);
	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30, 3);
	text = "power";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 5);

	text = "power";

	/* create filename with correct extension */
	plot(bitMap, text, "png");

	deleteBitMap(bitMap);

}


void Transient::plotPowerAss(int num_iter, double time){
	log_printf(NORMAL, "plotting Power ass...");

	int _bit_length_x = _plotter->getBitLengthX();
	int _bit_length_y = _plotter->getBitLengthY();

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _bit_length_x;
	bitMap->geom_y = _bit_length_y;
	bitMap->color_type = BLACKWHITE;

	std::stringstream text_stream;
	std::string text;
	std::string text2;

	/* draw and label axes */
	drawLine(bitMap, _bit_length_x / 10, 10, _bit_length_x / 10, _bit_length_y - 10);
	drawLine(bitMap, 10, 9 * _bit_length_y / 10, _bit_length_x - 10, 9 * _bit_length_y / 10);
	text = "power";
	drawText(bitMap, text, 30, _bit_length_y / 2 - 30);
	text = "time";
	drawText(bitMap, text, _bit_length_x / 2, _bit_length_y - 30);
	text_stream << time;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, 9.5 * _bit_length_x / 10, 9 * _bit_length_y / 10 + 20);


	/* create x axis scale */
	double power_max = 0, power_min = 2;
	for (int i = 1; i < num_iter; i++){
		if (_power_ass[i] > power_max){
			power_max = _power_ass[i];
		}
		if (_power_ass[i] > power_max){
			power_max = _power_ass[i];
		}
		if (_power_ass[i] < power_min){
			power_min = _power_ass[i];
		}
		if (_power_ass[i] < power_min){
			power_min = _power_ass[i];
		}
	}

	text_stream << power_min;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, 9 * _bit_length_y / 10 - 10);
	text_stream << power_max;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 10);
	double power_avg = (power_max - power_min)/2 + power_min;
	double power_range = power_max - power_min;
	text_stream << power_avg;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 2);

	text = "blue";
	text2 = "white";

	/* draw CMFD keff */
	for (int i = 1; i < num_iter; i++){
		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 2 - (_power_ass[i] - power_avg) / (power_range / 2) * (4 * _bit_length_y / 10), 3);
	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30, 3);
	text = "power";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 5);

	text = "power_ass";

	/* create filename with correct extension */
	plot(bitMap, text, "png");

	deleteBitMap(bitMap);

}


void Transient::plotTempMax(int num_iter, double time){
	log_printf(NORMAL, "plotting Temp max...");

	int _bit_length_x = _plotter->getBitLengthX();
	int _bit_length_y = _plotter->getBitLengthY();

	/* set up bitMap */
	BitMap<int>* bitMap = new BitMap<int>;
	bitMap->pixel_x = _bit_length_x;
	bitMap->pixel_y = _bit_length_y;
	initialize(bitMap);
	bitMap->geom_x = _bit_length_x;
	bitMap->geom_y = _bit_length_y;
	bitMap->color_type = BLACKWHITE;

	std::stringstream text_stream;
	std::string text;
	std::string text2;

	/* draw and label axes */
	drawLine(bitMap, _bit_length_x / 10, 10, _bit_length_x / 10, _bit_length_y - 10);
	drawLine(bitMap, 10, 9 * _bit_length_y / 10, _bit_length_x - 10, 9 * _bit_length_y / 10);
	text = "temp";
	drawText(bitMap, text, 30, _bit_length_y / 2 - 30);
	text = "time";
	drawText(bitMap, text, _bit_length_x / 2, _bit_length_y - 30);
	text_stream << time;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, 9.5 * _bit_length_x / 10, 9 * _bit_length_y / 10 + 20);


	/* create x axis scale */
	double temp_max = 0, temp_min = 2;
	for (int i = 1; i < num_iter; i++){
		if (_temp_max[i] > temp_max){
			temp_max = _temp_max[i];
		}
		if (_temp_max[i] > temp_max){
			temp_max = _temp_max[i];
		}
		if (_temp_max[i] < temp_min){
			temp_min = _temp_max[i];
		}
		if (_temp_max[i] < temp_min){
			temp_min = _temp_max[i];
		}
	}

	text_stream << temp_min;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, 9 * _bit_length_y / 10 - 10);
	text_stream << temp_max;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 10);
	double temp_avg = (temp_max - temp_min)/2 + temp_min;
	double temp_range = temp_max - temp_min;
	text_stream << temp_avg;
	text = text_stream.str();
	text_stream.str("");
	drawText(bitMap, text, _bit_length_x / 10 - 50, _bit_length_y / 2);

	text = "blue";
	text2 = "white";

	/* draw CMFD keff */
	for (int i = 1; i < num_iter; i++){
		drawPoint(bitMap, text, text2, 1, _bit_length_x / 10 + i * (8.5 * _bit_length_x / 10 / num_iter) , _bit_length_y / 2 - (_temp_max[i] - temp_avg) / (temp_range / 2) * (4 * _bit_length_y / 10), 3);
	}

	drawPoint(bitMap, text, text2, 1, 8 * _bit_length_x / 10, _bit_length_y / 30, 3);
	text = "Temp";
	drawText(bitMap, text, 8 * _bit_length_x / 10 + 20, _bit_length_y / 30 + 5);

	text = "temp_max";

	/* create filename with correct extension */
	plot(bitMap, text, "png");

	deleteBitMap(bitMap);

}

