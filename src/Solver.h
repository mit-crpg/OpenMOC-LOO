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
    int _cw;
    int _ch;
    int _num_crn;
    int _boundary_iteration;
    int _num_azim;
    int _num_FSRs;
    int* _num_tracks;
    double _total_vol;
    double _damp_factor;
    double _track_spacing;
    double _k_eff;
    double _cmfd_k;
    double _loo_k;
    double _k_half;
    double *_FSRs_to_fluxes[NUM_ENERGY_GROUPS + 1];
    double *_FSRs_to_powers;
    double *_FSRs_to_pin_powers;
    double *_FSRs_to_fission_source;
    double *_FSRs_to_scatter_source;
    double *_FSRs_to_absorption[NUM_ENERGY_GROUPS + 1];
    double *_FSRs_to_pin_absorption[NUM_ENERGY_GROUPS + 1];

    std::queue<double> _old_k_effs;
    std::queue<double> _old_eps_2;
    std::queue<double> _delta_phi;
    Plotter* _plotter;
    float* _pix_map_total_flux;
    Cmfd* _cmfd;
    std::string _geometry_file;

#if !STORE_PREFACTORS
    double* _pre_factor_array;
    int _pre_factor_array_size;
    int _pre_factor_max_index;
    double _pre_factor_spacing;
#endif
    double _l2_norm_conv_thresh;
    double _moc_conv_thresh;
    double computePreFactor(segment* seg, int energyg, int angle);
    bool _compute_powers;
    bool _run_cmfd;
    bool _run_loo;
    bool _run_loo1;
    bool _run_loo2;
    bool _diffusion;
    bool _acc_after_MOC_converge;
    bool _update_keff;
    bool _update_boundary;
    bool _plot_loo;
    void precomputeFactors();
    void initializeFSRs();
public:
    Solver(Geometry* geom, TrackGenerator* track_generator, 
           Plotter* plotter, Cmfd* cmfd, Options* opts);
    virtual ~Solver();

    /* setters */
    void setK(double k);
    void runCmfd();
    void runLoo1();
    void runLoo2();

    /* initialization */
    void initializeTrackFluxes(double flux);
    void oneFSRFluxOldSource();
    void zeroFSRFluxes();
    void zeroMeshCells();
    void zeroLeakage();
    void initializeWeights();

    /* checking and debugging */
    void checkTrackSpacing();
    void checkBoundary();
    void checkNeutronBalance();
    void checkNeutronBalanceWithDs();

    /* main routines */
    double kernel(int max_iterations);
    void MOCsweep(int max_iterations, int moc_iter);
    double runLoo(int moc_iter);
    double runCmfd(int moc_iter);
    void tallyLooWeight(Track *t, segment *seg, MeshSurface **surf, int dir);
    void tallyLooCurrent(Track *t, segment *seg, MeshSurface **surf, int dir);
    void tallyCmfdCurrent(Track *t, segment *seg, MeshSurface **surf, int dir);

    /* updates after transport sweep */
    void computeRatios();
    double computeKeff(int moc_iter);
    void normalizeFlux();
    void renormCurrents(Mesh* mesh, double keff);
    void updateSource();
    void updateFlux(int moc_iter);
    void updateBoundaryFluxByQuadrature();
    void storeFsrFluxPower();
    double computeFsrL2Norm(double *old_fsr_powers);
    double computeFsrLinf(double *old_fsr_powers);
    double computeSpectralRadius(double *old_fsr_powers); 
    void storeMOCBoundaryFlux();

    /* getters and setters */
    double** getFSRtoFluxMap();
    double getEps(Mesh* mesh, double keff, double renorm_factor);
    FlatSourceRegion* getFSRs();
    void setOldFSRFlux();

    /* printing and plotting */
    void printToScreen(int moc_iter);
    void printToLog(int moc_iter, double eps_inf, double eps_2, double rho);
    void printToMinimumLog(int moc_iter);
    void plotFluxes(int moc_iter);
    void plotPinPowers();
    void plotEverything(int moc_iter);
};

#endif /* SOLVER_H_ */
