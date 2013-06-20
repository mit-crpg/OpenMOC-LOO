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
#include <iostream>
#include <fstream>
#include "TrackGenerator.h"
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

#if USE_OPENMP
#include <omp.h>
#endif


/**
 * Solver types
 */
enum solveType {
	DIFFUSION,
	CMFD,
	LOO
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
	Mat _A;
	Mat _M;
	Vec _phi_new;
	Vec _source_old;
	double _damp;
	double _keff;
	double _l2_norm;
	double _l2_norm_conv_thresh;
	double _spacing;
	int _num_azim;
	bool _use_diffusion_correction;
	bool _run_loo;
	bool _run_cmfd;
	int _num_iter_to_conv;
public:
	Cmfd(Geometry* geom, Plotter* plotter, Mesh* mesh, 
		 bool runCmfd, bool runLoo, double damp,
		 bool useDiffusionCorrection, double l2_norm_conv_thresh,
		 TrackGenerator *track_generator);
	virtual ~Cmfd();
 	void computeDs();
	void computeDsBackup();
	void computeDsxDirection(double x, double y, int e, MeshCell *meshCell, 
							 double d, double f, double flux, int cell_width, 
							 double dt_weight);
 	void computeXS();
 	double computeDiffCorrect(double d, double h);
 	void updateMOCFlux(int iteration);
 	int constructAMPhi(Mat A, Mat B, Vec phi_old, solveType solveMethod);
 	double computeCMFDFluxPower(solveType solveMethod, int moc_iter);
 	Mat getA();
 	Mat getM();
 	Vec getPhiNew();
 	int createAMPhi(PetscInt size1, PetscInt size2, int cells);
 	double getKeff();
	void setOldFSRFlux();
	void setFSRs(FlatSourceRegion *fsrs);
	int fisSourceNorm(Vec snew, int iter, int num_cmfd_iteration);
	double getL2Norm();
	/* LOO */
	void storePreMOCMeshSource(FlatSourceRegion* fsrs);
	void computeQuadSrc();
	void computeQuadFlux();
	double computeLooFluxPower(solveType solveMethod, int moc_iter, double k);
	int getNumIterToConv();
};

#endif /* CMFD_H_ */
