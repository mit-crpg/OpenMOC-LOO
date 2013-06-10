/*
 * Solver.h
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <utility>
#include <math.h>
#include <unordered_map>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include "Geometry.h"
#include "Quadrature.h"
#include "Track.h"
#include "TrackGenerator.h"
#include "FlatSourceRegion.h"
#include "configurations.h"
#include "log.h"
#include "quickplot.h"
#include "Mesh.h"
#include "MeshCell.h"
#include "Material.h"
#include "Cmfd.h"
#include "Options.h"
#include "petsc.h"
#include <petscmat.h>

#if USE_OPENMP == true
	#include <omp.h>
#endif

class Solver {
private:
	Geometry* _geom;
	Quadrature* _quad;
	FlatSourceRegion* _flat_source_regions;
	Track** _tracks;
	int* _num_tracks;
	int _num_azim;
	int _num_FSRs;
	double *_FSRs_to_fluxes[NUM_ENERGY_GROUPS + 1];
	double *_FSRs_to_powers;
	double *_FSRs_to_pin_powers;
	double *_FSRs_to_fission_source;
	double *_FSRs_to_scatter_source;
	double *_FSRs_to_absorption[NUM_ENERGY_GROUPS + 1];
	double *_FSRs_to_pin_absorption[NUM_ENERGY_GROUPS + 1];
	double _k_eff;
	double _cmfd_k;
	double _loo_k;
	std::queue<double> _old_k_effs;
	Plotter* _plotter;
	float* _pix_map_total_flux;
	Cmfd* _cmfd;
#if !STORE_PREFACTORS
	double* _pre_factor_array;
	int _pre_factor_array_size;
	int _pre_factor_max_index;
	double _pre_factor_spacing;
	bool _update_flux;
	double _l2_norm_conv_thresh;
	double _moc_conv_thresh;
	bool _compute_powers;
	bool _run_cmfd;
	bool _run_loo;
	bool _diffusion;
	bool _loo_after_MOC_converge;
#endif
	void precomputeFactors();
	double computePreFactor(segment* seg, int energyg, int angle);
	void initializeFSRs();
public:
	Solver(Geometry* geom, TrackGenerator* track_generator, 
		   Plotter* plotter, Cmfd* cmfd, Options* opts);
/*		   bool _update_flux, double l2NormConvThresh, 
		   bool computePowers, bool runCmfd, bool runLoo, bool diffusion, 
		   double k_guess); */
	virtual ~Solver();
	void initializeTrackFluxes(double flux);
	void initializeSource();
	void normalizeFlux();
 	void updateSource();
	void oneFSRFluxes();
	void zeroFSRFluxes();
	void zeroMeshCells();
	void zeroLeakage();
	void computeRatios();
	void updateKeff(int iteration);
	double** getFSRtoFluxMap();
	void MOCsweep(int max_iterations);
	double kernel(int max_iterations);
	void plotFluxes(int iter_num);
	void checkTrackSpacing();
	void computeFsrPowers();
	void plotPinPowers();
 	void checkNeutBal(Mesh* mesh);
 	void renormCurrents(Mesh* mesh, double keff);
 	double getEps(Mesh* mesh, double keff, double renorm_factor);
 	FlatSourceRegion* getFSRs();
 	void setOldFSRFlux();
	void tallyLooCurrent(Track *track, segment *segment, 
						 MeshSurface **meshSurfaces, int dir);
	void tallyCmfdCurrent(Track *track, segment *segment, 
						 MeshSurface **meshSurfaces, int dir);
	double runLoo(int i);
	double runCmfd(int i);
	double computeL2Norm(double *old_fsr_powers);
};

#endif /* SOLVER_H_ */
