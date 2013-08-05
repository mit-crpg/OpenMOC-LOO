/*
 * Options.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#include "Options.h"

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)


/**
 * Options constructor
 * @param argc the number of command line arguments from console
 * @param argv a char array of command line arguments from console
 */
Options::Options(int argc, char **argv) 
{
    /* convergence options */
    _l2_norm_conv_thresh = 1e-10; /* convergence on acceleration iteration */
    _moc_conv_thresh = 1e-8;      /* convergence on MOC sweeps */

    /* most important acceleration options */
    _acc_after_MOC_converge = false;
    _cmfd = false; 					
    _loo = false;
    _loo1 = false;
    _loo2 = false;
    _damp_factor = 1.0;
    _boundary_iteration = 0;
    _diffusion = true;			/* run diffusion for 1st iter */
    _update_boundary = true;            /* update boundary angular flux */

    if (std::string(getenv("PWD")).find("Release") != std::string::npos)
        _relative_path = "../";
    else
        _relative_path = "";

    /* Default geometry input file */
    _geometry_file = _relative_path + 
        "xml-sample/Cmfd/geometry_homo.xml";
        //"xml-sample/Cmfd/geometry_8x8_leakage3_2.xml"; 
    _material_file = _relative_path + "xml-sample/Cmfd/material_simple.xml";
	
    _track_spacing = 0.5;		/* Default track spacing in C4: 0.5cm */
    _num_azim = 64;			/* Default # azimuthal angle in C4: 64*/
    //_track_spacing = 0.8909545443;
    //_num_azim = 4;

    /* MOC options */
    _verbosity = "NORMAL";		/* Default logging level */
    _dump_geometry = false;		/* Default: not dump geometry */
    _extension = "png";			/* Default: plot in png */
    _compute_pin_powers = false;	/* Default: not compute pin powers */
    _compress_cross_sections = false;	/* Default: not compress xs */

    /* other acceleration options */
    _cmfd_level = 0;			/* Default cmfd level is 0: pin-wise */
    _multigroup = true;			/* sets CMFD to one group structure */
    _k_guess = 1.0;
    _update_keff = false;  		/* Default: not use CMFD's k */
    _print_matrices = false;		/* Default will not print matrices */
    _diffusion_correction = false; 

    /* plotting options */
    _bit_dimension = 1000;		/* y dimension of plots */
    _plot_specs = false;            	/* plot mat, cell, FSR, track, seg */
    _plot_quad_flux = false;            /* plot quad flux, net current, xs */
    _plot_fluxes = false;		/* plot colors, not values*/
    _plot_current = false;		/* plot cmfd currents */
    _plot_diffusion = false;		/* FIXME: does nothing now */
    _plot_keff = false;			/* Default will not plot keff */
    _plot_prolongation = false;

    /* All extra options get placed into this array, which can be used
     * to call sub-initializers (petsc, for instance) */
    extra_argc = 0;
    extra_argv = (char **)malloc(sizeof(*this->extra_argv) * argc);
    for (int i = 0 ; i < argc; i++)
        extra_argv[i] = NULL;

    for (int i = 0; i < argc; i++) {
        if (i > 0) {
            if (LAST("--geometryfile") || LAST("-g")) {
                _geometry_file = argv[i];
            }
            else if (LAST("--materialfile") || LAST("-m")) {
                _material_file = argv[i];
            }
            else if (LAST("--trackspacing") || LAST("-ts"))
                _track_spacing = atof(argv[i]);
            else if (LAST("--numazimuthal") || LAST("-na"))
                _num_azim = atoi(argv[i]);
            else if (LAST("--boundaryiteration") || LAST("-bi"))
                _boundary_iteration = atoi(argv[i]);
            else if (LAST("--dampfactor") || LAST("-df"))
                _damp_factor = atof(argv[i]);
            else if (LAST("--bitdimension") || LAST("-bd"))
                _bit_dimension = atoi(argv[i]);
            else if (LAST("--verbosity") || LAST("-v"))
                _verbosity = strdup(argv[i]);
            else if (LAST("--kguess") || LAST("-k"))
                _k_guess = atof(argv[i]);
            else if (strcmp(argv[i], "-dg") == 0 ||
                     strcmp(argv[i], "--dumpgeometry") == 0)
                _dump_geometry = true;
            else if (LAST("--extension") || LAST("-ex"))
                _extension = argv[i];
            else if (strcmp(argv[i], "-noconv") == 0)
                _acc_after_MOC_converge = false;
            else if (strcmp(argv[i], "-debug") == 0)
                _acc_after_MOC_converge = true;
            else if (strcmp(argv[i], "-ps") == 0 ||
                     strcmp(argv[i], "--plotspecs") == 0)
                _plot_specs = true;
            else if (strcmp(argv[i], "-pf") == 0 ||
                     strcmp(argv[i], "--plotfluxes") == 0)
                _plot_fluxes = true;
            else if (strcmp(argv[i], "-pp") == 0 ||
                     strcmp(argv[i], "--plotprolongation") == 0)
                _plot_prolongation = true;
            else if (strcmp(argv[i], "-cp") == 0 ||
                     strcmp(argv[i], "--computepowers") == 0)
                _compute_pin_powers = true;
            else if (strcmp(argv[i], "-cxs") == 0 ||
                     strcmp(argv[i], "--compressxs") == 0)
                _compress_cross_sections = true;
            else if (strcmp(argv[i], "-uk") == 0 ||
                     strcmp(argv[i], "--updatekeff") == 0)
                _update_keff = true;
            else if (strcmp(argv[i], "-ub") == 0 ||
                     strcmp(argv[i], "--updateboundary") == 0)
                _update_boundary = true;
            else if (strcmp(argv[i], "-nub") == 0 ||
                     strcmp(argv[i], "--noupdateboundary") == 0)
                _update_boundary = false;
            else if (strcmp(argv[i], "-nc") == 0 ||
                     strcmp(argv[i], "--nocmfd") == 0)
                _cmfd = false;
            else if (strcmp(argv[i], "-wc") == 0 ||
                     strcmp(argv[i], "--withcmfd") == 0)
            {
                _cmfd = true;
                _loo = false;
                _damp_factor = 0.66;
            }
            else if (strcmp(argv[i], "-nl") == 0 ||
                     strcmp(argv[i], "--noloo") == 0)
                _loo = false;
            else if (strcmp(argv[i], "-wl") == 0 ||
                     strcmp(argv[i], "--withloo") == 0)
            {
                _cmfd = false;
                _loo = true;
                _loo1 = false;
                _loo2 = true;
                _damp_factor = 1.0;
            }
            else if (strcmp(argv[i], "-wl1") == 0 ||
                     strcmp(argv[i], "--withloo1") == 0)
            {
                _cmfd = false;
                _loo = true;
                _loo1 = true;
                _loo2 = false;
                _damp_factor = 1.0;
            }
            else if (strcmp(argv[i], "-wl2") == 0 ||
                     strcmp(argv[i], "--withloo2") == 0)
            {
                _cmfd = false;
                _loo = true;
                _loo1 = false;
                _loo2 = true;
                _damp_factor = 1.0;
            }
            else if (strcmp(argv[i], "-pc") == 0 ||
                     strcmp(argv[i], "--plotcurrent") == 0)
                _plot_current = true;
            else if (strcmp(argv[i], "-pk") == 0 ||
                     strcmp(argv[i], "--plotkeff") == 0)
                _plot_keff = true;
            else if (strcmp(argv[i], "-diff") == 0 ||
                     strcmp(argv[i], "--diffusion") == 0)
                _diffusion = true;
            else if (strcmp(argv[i], "-ndiff") == 0 ||
                     strcmp(argv[i], "--nodiffusion") == 0)
                _diffusion = false;
            else if (strcmp(argv[i], "-pd") == 0 ||
                     strcmp(argv[i], "--plotdiffusion") == 0)
                _plot_diffusion = true;
            else if (LAST("--fluxconv") || LAST("-fc"))
            {
                _moc_conv_thresh = atof(argv[i]);
                _l2_norm_conv_thresh = _moc_conv_thresh * 1e-2;
            }
            else if (LAST("--l2normconv") || LAST("-lc"))
                _l2_norm_conv_thresh = atof(argv[i]);
            else if (strcmp(argv[i], "-mg") == 0 ||
                     strcmp(argv[i], "--multigroup") == 0)
                _multigroup = true;
            else if (strcmp(argv[i], "-sg") == 0 ||
                     strcmp(argv[i], "--singlegroup") == 0)
                _multigroup = false;			
            else if (strcmp(argv[i], "-pm") == 0 ||
                     strcmp(argv[i], "--printmatrices") == 0)
                _print_matrices = true;
            else if (LAST("--cmfdlevel") || LAST("-cl"))
                _cmfd_level = atoi(argv[i]);
            else if (strcmp(argv[i], "-dc")==0 || 
                     strcmp(argv[i], "--diffusioncorrection") == 0)
                _diffusion_correction = true;
            else
                this->extra_argv[this->extra_argc++] = strdup(argv[i]);
        }
    }

    /* If debug mode is on for CMFD, damping factor has to be 1.0 because 
     * otherwise D tilde is picking up 0.0 as old values */
    if ((_acc_after_MOC_converge) && (_cmfd))
        _damp_factor = 1.0;

}

Options::~Options(void) { }

/**
 * Returns a character array with the path to the geometry input file. By 
 * default this will return the path to /xml-sample/1/geometry.xml if not set 
 * at runtime from the console
 * @return path to the geometry input file
 */
const char *Options::getGeometryFile() const {
    return (const char*)_geometry_file.c_str();
}

/**
 * Returns a character array with the path to the material input file. By 
 * default this will return the path to /xml-sample/1/material.xml if not set 
 * at runtime from the console
 * @return path to the geometry input file
 */
const char *Options::getMaterialFile() const {
    return (const char*)_material_file.c_str();
}

/**
 * Returns a boolean representing whether or not to dump the geometry to the
 * console. If true, the geometry will be printed out after parsing is complete
 * @return whether or not to dump the geometry to the console
 */
bool Options::dumpGeometry(){
    return _dump_geometry;
}

/**
 * Returns the number of azimuthal angles. 
 * @return the number of azimuthal angles
 */
double Options::getNumAzim(){
    return _num_azim;
}

/**
 * Returns the y dimension of plots. By default this will return 1000 bits
 * (or pixels) if not set at runtime from the console
 * @return the y dimension of plots.
 */
int Options::getBitDimension(){
    return _bit_dimension;
}



/**
 * Returns the track spacing. By default this will return 0.05 if not set 
 * at runtime from the console
 * @return the track spacing
 */
double Options::getTrackSpacing(){
    return _track_spacing;
}

double Options::getKGuess(){
    return _k_guess;
}

/**
 * Returns the verbosity logging level. By default this will return NORMAL if 
 * not set at runtime from the console
 * @return the verbosity
 */
char* Options::getVerbosity(){
    return (char*)_verbosity.c_str();
}

/**
 * Returns the image files extension. By default this will return .png if not 
 * set at runtime from the console
 * @return the image files extension
 */
std::string Options::getExtension(){
    return _extension;
}

std::string Options::getGeometryFile(){
    return _geometry_file;
}

/**
 * Returns a boolean representing whether or not to plot the specs.
 *  If true, the specs will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotSpecs(){
    return _plot_specs;
}

/**
 * Returns a boolean representing whether or not to plot the cells.
 *  If true, the cells will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotFluxes(){
    return _plot_fluxes;
}


/**
 * Returns a boolean representing whether or not to compute the powers
 * in each pin. If true, txt files with the pin powers will be created
 * in a new directory called "PinPowers"
 * @return whether or not to compute the pin powers
 */
bool Options::computePinPowers(){
    return _compute_pin_powers;
}


/**
 * Returns a boolean representing whether or not to compress the
 * cross-sections for each material. If true, the starting and ending
 * index for the non-zero cross-section values will be computed to
 * help speed up fixed source iteraiton. Note: this will only speed
 * up fixed source iteration for materials with many zeroes in their
 * cross-section values
 * @return whether or not to compute the pin powers
 */
bool Options::compressCrossSections(){
    return _compress_cross_sections;
}

/**
 * Returns a boolean representing whether or not to perform CMFD acceleration
 * @return whether or not to perform CMFD acceleration
 */
bool Options::cmfd(){
    return _cmfd;
}

/**
 * Returns a boolean representing whether or not to plot the cmfd fluxes
 * at each step. If true, the net current will be plotted in a file of 
 * _extension type
 * @return whether or not to plot net current
 */
bool Options::plotCurrent(){
    return _plot_current;
}

/**
 * Returns a boolean representing whether or not to plot the loo surface
 * averaged quadrature fluxes at each step. If true, the quadrature fluxes 
 * will be plotted in a file of _extension type
 * @return whether or not to plot quadrature flux
 */
bool Options::plotQuadFlux(){
    return _plot_quad_flux;
}

/**
 * Returns a boolean representing whether or not to plot the diffusion flux.
 *  If true, the net current will be plotted in a file of _extension type
 * @return whether or not to plot net current
 */
bool Options::plotDiffusion(){
    return _plot_diffusion;
}

/**
 * Returns a boolean representing whether or not to plot keff.
 *  If true, the net current will be plotted in a file of _extension type
 * @return whether or not to plot net current
 */
bool Options::plotKeff(){
    return _plot_keff;
}

/**
 * Returns a boolean representing whether or not to use CMFD to update flux
 * @return whether or not to use CMFD to update flux
 */
bool Options::updateKeff(){
    return _update_keff;
}

/**
 * Returns the fission source convergence threshold for acceleration iteration
 * @return fission source convergence threshold innter iteration
 */
double Options::getL2NormConvThresh() {
    return _l2_norm_conv_thresh;
}

/**
 * Returns the fission source convergence threshold for MOC outter iteration
 * @return fission source convergence threshold for outter iteration
 */
double Options::getMOCConvThresh() {
    return _moc_conv_thresh;
}


/**
 * Returns bool telling us cmfd group structure
 * @return bool telling us cmfd group structure
 */
bool Options::getGroupStructure() {
    return _multigroup;
}


/**
 * Returns bool telling us cmfd group structure
 * @return bool telling us cmfd group structure
 */
bool Options::getPrintMatrices() {
    return _print_matrices;
}

int Options::getCmfdLevel(){
    return _cmfd_level;
}

bool Options::getCmfd(){
    return _cmfd;
}

bool Options::getLoo(){
    return _loo;
}

bool Options::getLoo1(){
    return _loo1;
}

bool Options::getLoo2(){
    return _loo2;
}

bool Options::getDiffusion(){
    return _diffusion;
}

bool Options::getDiffusionCorrection(){
    return _diffusion_correction;
}

bool Options::getAccAfterMOCConverge(){
    return _acc_after_MOC_converge;
}

bool Options::plotProlongation(){
    return _plot_prolongation;
}

double Options::getDampFactor()
{
    return _damp_factor;
}

int Options::getBoundaryIteration()
{
    return _boundary_iteration;
}

bool Options::getUpdateBoundary()
{
    return _update_boundary;
}
