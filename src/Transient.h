/*
 * Transient.h
 *
 *  Created on: Dec 15, 2012
 *      Author: samuelshaner
 */

#ifndef TRANSIENT_H_
#define TRANSIENT_H_

#include <math.h>
#include <string>
#include <sstream>
#include "Geometry.h"
#include "Quadrature.h"
#include "FlatSourceRegion.h"
#include "configurations.h"
#include "log.h"
#include "quickplot.h"
#include <armadillo>
#include "Cmfd.h"
#include "Mesh.h"
#include "MeshCell.h"
#include "MeshSurface.h"
#include "Solver.h"
#include "Plotter.h"
#include "quickplot.h"

class Transient {
private:
	Geometry* _geom;
	Cmfd* _cmfd;
	Mesh* _mesh;
	Solver* _solver;
	Plotter* _plotter;
	double _precursors[2];
	double _amp[2];
	double _beta[2];
	double _velocities[2];
	double _lambda[2];
	double _sigma_A_1;
	double _sigma_S_4;
	double _alpha;
	double _gamma;
	double _kappa;
	double _nu;
	double _buckle;
	double _temp_0;
	double _t_max;
	double _dt_outer;
	double _dt_inner;
	int _t_steps;
	double _power_init;
	double _vol_core;
	double _keff_0;
	double *_temp;
	double *_power;
	double *_power_pke;
	double *_power_ass;
	double *_temp_max;
	double _power_factor;

public:
	Transient(Geometry* geom, Cmfd* cmfd, Mesh* mesh, Solver* solver, Plotter* plotter, double tMax, double deltaT, double keff_0);
	virtual ~Transient();
	void solve();
	int setNMatrix(Mat A, double keff_0);
	double computePower(Mesh* mesh, arma::colvec B);
	double computePowerAss(Mesh* mesh, arma::colvec B);
	double computeTemp(Mesh* mesh, arma::colvec B);
	double computeTempMax(Mesh* mesh, arma::colvec B);
	void computeInitialPrecursorConc(Mesh* mesh);
	arma::mat constructA(Mesh* mesh, arma::mat A, solveType solveMethod);
	arma::mat constructM(Mesh* mesh, arma::mat M);
	arma::mat constructN(Mesh* mesh, arma::mat N, arma::mat A, arma::mat M, arma::colvec shape);
	arma::colvec constructShape(arma::colvec shape, Mesh* mesh);
	arma::colvec constructB(arma::colvec B, arma::mat N, arma::mat M, arma::colvec shape);
	int normalizeN(Mat N, Vec shape);
	Mesh* getMesh();
	void copyMesh(Mesh* mesh1, Mesh* mesh2);
	void moveBlade(Mesh* mesh, double time);
	void tempFeedback(Mesh* mesh, arma::colvec B, double time_new);
	void plotTemp(int num_iter, double time);
	void plotPower(int num_iter, double time);
	void plotPowerAss(int num_iter, double time);
	void plotTempMax(int num_iter, double time);
	void plotTempAss(Mesh* mesh);
};


#endif /* TRANSIENT_H_ */
