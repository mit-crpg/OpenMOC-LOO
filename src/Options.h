/*
 * Options.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string.h>
#include <stdlib.h>
#include "log.h"

class Options {
private:
    std::string _relative_path;
    std::string _geometry_file;
    std::string _material_file;
    std::string _track_input_file;
    std::string _extension;
    std::string _verbosity;
    int _num_azim;
    int _bit_dimension;
    int _cmfd_level;
    int _boundary_iteration;
    double _track_spacing;
    double _l2_norm_conv_thresh;
    double _moc_conv_thresh;
    double _k_guess;
    double _damp_factor;
    bool _dump_geometry;
    bool _plot_specs;
    bool _plot_fluxes;
    bool _compute_pin_powers;
    bool _compress_cross_sections;
    bool _run_all;
    bool _cmfd;
    bool _loo;
    bool _loo1; //update by psi
    bool _loo2; //update by phi
    bool _acc_after_MOC_converge;
    bool _plot_quad_flux;
    bool _plot_current;
    bool _plot_diffusion;
    bool _plot_keff;
    bool _plot_prolongation;
    bool _update_keff;
    bool _multigroup;
    bool _print_matrices;
    bool _diffusion;
    bool _diffusion_correction;
    bool _update_boundary;
    bool _reflect_outgoing;
public:
    Options(int argc, char **argv);
    ~Options(void);
    std::string getGeometryFile();
    std::string getExtension();
    char **extra_argv;
    char* getVerbosity();
    const char *getGeometryFile() const;
    const char *getMaterialFile() const;

    int extra_argc;
    int getBitDimension();
    int getCmfdLevel();
    int getBoundaryIteration();
    double getTrackSpacing();
    double getNumAzim();
    double getL2NormConvThresh();
    double getMOCConvThresh();
    double getKGuess();
    double getDampFactor();
    bool dumpGeometry();
    bool plotSpecs();
    bool plotFluxes();
    bool computePinPowers();
    bool compressCrossSections();
    bool cmfd();
    bool plotQuadFlux();
    bool plotCurrent();
    bool plotDiffusion();
    bool plotKeff();
    bool plotProlongation();
    bool updateKeff();
    bool getGroupStructure();
    bool getPrintMatrices();
    bool getRunAll();
    bool getCmfd();
    bool getLoo();
    bool getLoo1();
    bool getLoo2();
    bool getDiffusion();
    bool getDiffusionCorrection();
    bool getAccAfterMOCConverge();
    bool getUpdateBoundary();
    bool getReflectOutgoing();
};

#endif
