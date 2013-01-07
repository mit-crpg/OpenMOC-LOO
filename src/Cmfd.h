/*
 * Cmfd.h
 *
 *  Created on: September 13, 2012
 *      Author: Sam Shaner
 *				MIT, Course 22
 *              shaner@mit.edu
 */

#ifndef CMFD_H_
#define CMFD_H_

#include <utility>
#include <math.h>
#include <unordered_map>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include "Geometry.h"
#include "Quadrature.h"
#include "FlatSourceRegion.h"
#include "configurations.h"
#include "log.h"
#include "quickplot.h"
#include "Mesh.h"
#include "MeshCell.h"
#include "MeshSurface.h"
#include "Material.h"
#include "Surface.h"
#include "petsc.h"
#include <petscmat.h>

/**
 * Surface types
 */
enum solveType {
	DIFFUSION,
	CMFD
};

#include "Plotter.h"

class Cmfd {
private:
	Geometry* _geom;
	Quadrature* _quad;
	Mesh* _mesh;
	FlatSourceRegion* _flat_source_regions;
	double _k_eff;
	Plotter* _plotter;
	bool _update_flux;
	Mat _A;
	Mat _M;
	Vec _phi_new;
	double _keff;

public:
	Cmfd(Geometry* geom, Plotter* plotter, Mesh* mesh, bool updateFlux, bool runCmfd);
	virtual ~Cmfd();
 	void computeDs();
 	void computeXS(FlatSourceRegion* fsrs);
 	double computeDiffCorrect(double d, double h);
 	void updateMOCFlux();
 	int constructAMPhi(Mat A, Mat B, Vec phi_old, solveType solveMethod);
 	double computeCMFDFluxPower(solveType solveMethod, int moc_iter);
 	Mat getA();
 	Mat getM();
 	Vec getPhiNew();
 	int createAMPhi(PetscInt size1, PetscInt size2, int cells);
 	double getKeff();


};

#endif /* CMFD_H_ */
